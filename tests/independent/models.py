"""Model bank of the independent review: raw benchfea input dictionaries (first order).

``LINEAR[name] = builder`` returns ``(raw, combination id)``.  Sections are of type ``properties`` so
that no section routine of benchfea is involved.  The models are deliberately irregular (skew
members, rolls, unequal spans, mixed unit systems, both vertical axes).
"""
from __future__ import annotations

from harness import base, case, member, nodes

STEEL = {("kN", "m"): {"E": 2.0e8, "G": 7.7e7, "unitWeight": 77.0},
         ("kN", "mm"): {"E": 200.0, "G": 77.0, "unitWeight": 7.7e-8},
         ("N", "mm"): {"E": 2.0e5, "G": 7.7e4, "unitWeight": 7.7e-5},
         ("kip", "in"): {"E": 29000.0, "G": 11200.0, "unitWeight": 2.84e-4},
         ("kip", "ft"): {"E": 4176000.0, "G": 1612800.0, "unitWeight": 0.49},
         ("lbf", "ft"): {"E": 4.176e9, "G": 1.6128e9, "unitWeight": 490.0}}
# A, Iy, Iz, J (+ shear areas) in m, scaled by the length unit
_SEC_M = {"A": 0.012, "Iy": 4.0e-5, "Iz": 1.5e-4, "J": 3.0e-6}
_LEN = {"m": 1.0, "mm": 1000.0, "in": 1 / 0.0254, "ft": 1 / 0.3048}


def sec(units, shear=False, scale=1.0, **over):
    f = _LEN[units[1]]
    s = {"type": "properties", "A": _SEC_M["A"] * scale * f ** 2, "Iy": _SEC_M["Iy"] * scale ** 2 * f ** 4,
         "Iz": _SEC_M["Iz"] * scale ** 2 * f ** 4, "J": _SEC_M["J"] * scale ** 2 * f ** 4}
    if shear:
        s["Asy"], s["Asz"] = 0.4 * s["A"], 0.6 * s["A"]
    s.update(over)
    return s


def start(units=("kN", "m"), vertical="Z", shear=False, **analysis):
    d = base(units, vertical, shearDeformation=shear, **analysis)
    d["materials"]["M"] = dict(STEEL[units])
    d["sections"]["S"] = sec(units, shear)
    d["sections"]["B"] = sec(units, shear, scale=2.0)
    return d


LINEAR = {}


def model(f):
    LINEAR[f.__name__] = f
    return f


@model
def L01_skew_rolls_everything():
    d = start(shear=True)
    nodes(d, A=(0, 0, 0), B=(3, 1, 4), C=(7, 2, 4.5), D=(7, 5, 0))
    d["supports"] = {"A": {"type": "fixed"}, "D": {"type": "pinned", "dofs": {"rx": {"spring": 5000}}}}
    member(d, "m1", "A", "B", orientation={"roll": 25})
    member(d, "m2", "B", "C", releases={"i": {"rz": {"spring": 3000}}, "j": {"ry": True}},
           offsets={"i": [0.2, 0.1, -0.05]})
    member(d, "m3", "C", "D", orientation={"roll": -40})
    case(d, "LC", {"type": "selfWeight"}, {"type": "nodal", "node": "B", "Fx": 10, "Fy": -5, "Mz": 3},
         {"type": "memberDistributed", "member": "m2", "direction": "Z", "w1": -4, "w2": -9, "start": 0.2,
          "end": 0.8, "projected": True},
         {"type": "memberPoint", "member": "m3", "direction": "Y", "value": 7, "position": 0.3},
         {"type": "memberPoint", "member": "m1", "direction": "x", "kind": "moment", "value": 2, "position": 0.6})
    return d, "LC"


@model
def L02_space_frame_y_up_kip_in():
    d = start(("kip", "in"), "Y")
    nodes(d, a=(0, 0, 0), b=(240, 0, 30), c=(250, 0, 200), e=(10, 0, 180),
          A=(0, 144, 0), B=(240, 150, 30), C=(250, 160, 200), E=(10, 140, 180))
    for n in "abce":
        d["supports"][n] = {"type": "fixed" if n in "ac" else "pinned"}
    for k, (lo, up) in enumerate(zip("abce", "ABCE")):
        member(d, f"col{k}", lo, up, orientation={"roll": 15 * k})
    for k, (i, j) in enumerate((("A", "B"), ("B", "C"), ("C", "E"), ("E", "A"))):
        member(d, f"bm{k}", i, j, sec="B", orientation={"roll": 10 * k - 12})
    member(d, "diag", "A", "C", orientation={"referenceVector": [0.2, 1, 0.1]})
    case(d, "G", *[{"type": "memberDistributed", "member": f"bm{k}", "direction": "Y", "w1": -0.08 - 0.01 * k}
                   for k in range(4)],
         {"type": "nodal", "node": "B", "Fx": 12, "Fz": -7, "My": 150},
         {"type": "memberPoint", "member": "diag", "direction": "Y", "value": -9, "position": 100,
          "positionType": "absolute"},
         {"type": "memberDistributed", "member": "col1", "direction": "z", "w1": 0.03, "start": 0.25, "end": 0.75})
    return d, "G"


@model
def L03_partial_releases_axis_aligned_kn_mm():
    d = start(("kN", "mm"))
    nodes(d, n1=(0, 0, 0), n2=(0, 0, 3500), n3=(6000, 0, 3500), n4=(6000, 0, 0), n5=(6000, 4000, 3500),
          n6=(6000, 4000, 0))
    for n in ("n1", "n4", "n6"):
        d["supports"][n] = {"type": "fixed"}
    member(d, "c1", "n1", "n2")
    member(d, "c2", "n4", "n3", releases={"j": {"ux": {"spring": 250.0}}})
    member(d, "c3", "n6", "n5", releases={"i": {"rx": {"spring": 4.0e5}, "ry": {"spring": 2.0e6}}})
    member(d, "b1", "n2", "n3", sec="B", releases={"i": {"rz": {"spring": 3.0e6}}, "j": {"rz": {"spring": 8.0e5},
                                                                                         "uy": {"spring": 40.0}}})
    member(d, "b2", "n3", "n5", sec="B", releases={"i": {"ry": {"spring": 5.0e5}, "uz": {"spring": 15.0}},
                                                   "j": {"rz": True}})
    case(d, "Q", {"type": "memberDistributed", "member": "b1", "direction": "Z", "w1": -0.02},
         {"type": "memberDistributed", "member": "b2", "direction": "y", "w1": -0.015, "start": 0.3},
         {"type": "memberDistributed", "member": "b2", "direction": "X", "w1": 0.004},
         {"type": "nodal", "node": "n2", "Fx": 30, "Fy": 12},
         {"type": "memberPoint", "member": "c3", "direction": "X", "value": 8, "position": 0.4})
    return d, "Q"


@model
def L04_full_releases_of_every_kind():
    d = start()
    nodes(d, n1=(0, 0, 0), n2=(4, 0, 0), n3=(8, 0, 0), n4=(8, 5, 0), n5=(8, 5, 3), n6=(4, 0, 3), n7=(12, 5, 3))
    d["supports"] = {"n1": {"type": "fixed"}, "n4": {"type": "fixed"}, "n6": {"type": "fixed"},
                     "n7": {"type": "fixed"}}
    member(d, "a", "n1", "n2", releases={"j": {"ux": True}})
    member(d, "b", "n2", "n3", releases={"i": {"rz": True}, "j": {"uz": True}})
    member(d, "c", "n3", "n4", releases={"i": {"rx": True, "ry": True}})
    member(d, "e", "n2", "n6", releases={"j": {"rz": True, "ry": True}})
    member(d, "f", "n4", "n5", releases={"j": {"uz": True}})
    member(d, "g", "n5", "n7", releases={"i": {"rx": True}})
    case(d, "P", {"type": "nodal", "node": "n2", "Fx": 20, "Fy": 6, "Fz": -15, "Mx": 2},
         {"type": "nodal", "node": "n3", "Fz": -25, "Fx": 5, "My": 4},
         {"type": "nodal", "node": "n5", "Fx": 3, "Fy": 4, "Fz": -6, "Mx": 1.5},
         {"type": "memberDistributed", "member": "b", "direction": "Z", "w1": -6},
         {"type": "memberDistributed", "member": "c", "direction": "Z", "w1": -3, "start": 0.2, "end": 0.7},
         {"type": "memberPoint", "member": "a", "direction": "y", "value": -9, "position": 0.35})
    return d, "P"


@model
def L05_offsets_off_axis_skew():
    d = start()
    nodes(d, A=(0, 0, 0), B=(0.5, 0.3, 4), C=(6, 1.5, 5), D=(6.5, 1, 0), F=(3, 6, 4.5), H=(3, 6.5, 0))
    for n in "ADH":
        d["supports"][n] = {"type": "fixed"}
    member(d, "c1", "A", "B", offsets={"j": [0, 0, -0.3]})
    member(d, "c2", "D", "C", orientation={"roll": 30}, offsets={"axes": "local", "j": [-0.35, 0.05, 0.02]})
    member(d, "c3", "H", "F")
    member(d, "b1", "B", "C", sec="B", orientation={"roll": 8}, offsets={"i": [0.25, 0.05, 0.1], "j": [-0.3, 0.0, -0.15]})
    member(d, "b2", "C", "F", sec="B", offsets={"axes": "local", "i": [0.2, 0.12, -0.08], "j": [-0.15, -0.1, 0.05]})
    member(d, "b3", "F", "B", sec="B", offsets={"i": [0, 0, 0.2], "j": [0.1, 0.2, 0.0]})
    case(d, "W", {"type": "memberDistributed", "member": "b1", "direction": "Z", "w1": -12},
         {"type": "memberDistributed", "member": "b2", "direction": "y", "w1": -7, "start": 0.1, "end": 0.9},
         {"type": "memberDistributed", "member": "b3", "direction": "X", "w1": 2.5},
         {"type": "memberPoint", "member": "b3", "direction": "Z", "value": -14, "position": 0.45},
         {"type": "nodal", "node": "C", "Fx": 9, "Fy": -11, "Mz": 6}, {"type": "selfWeight"})
    return d, "W"


@model
def L06_releases_together_with_offsets():
    d = start(("kip", "ft"), "Y")
    nodes(d, A=(0, 0, 0), B=(0, 12, 0), C=(20, 13, 2), D=(21, 0, 3), F=(20, 13, 18), H=(20, 0, 18))
    for n in "ADH":
        d["supports"][n] = {"type": "fixed"}
    member(d, "c1", "A", "B")
    member(d, "c2", "D", "C")
    member(d, "c3", "H", "F")
    member(d, "b1", "B", "C", sec="B", offsets={"i": [0.8, 0, 0.1], "j": [-0.7, 0.2, 0]},
           releases={"i": {"rz": True}, "j": {"rz": True, "ry": True}})
    member(d, "b2", "C", "F", sec="B", offsets={"axes": "local", "i": [0.9, 0.3, 0], "j": [-0.5, 0, 0.2]},
           releases={"i": {"rz": {"spring": 9000.0}, "ux": {"spring": 700.0}}, "j": {"rx": True, "rz": {"spring": 100.0}}})
    case(d, "W", {"type": "memberDistributed", "member": "b1", "direction": "Y", "w1": -1.2},
         {"type": "memberDistributed", "member": "b2", "direction": "Y", "w1": -0.5, "w2": -1.5},
         {"type": "memberPoint", "member": "b2", "direction": "x", "value": 6, "position": 0.6},
         {"type": "nodal", "node": "C", "Fx": 5, "Fz": 3})
    return d, "W"


@model
def L07_rigid_links_with_lever_arms():
    d = start()
    nodes(d, A=(0, 0, 0), B=(0, 0, 4), C=(5, 1, 4), D=(5, 1, 0), M=(2.5, 0.5, 4.6), S1=(2.0, 2.5, 5.5),
          S2=(3.4, -1.5, 3.8), T=(6.5, 3.0, 5.0))
    d["supports"] = {"A": {"type": "fixed"}, "D": {"type": "pinned"}}
    member(d, "c1", "A", "B")
    member(d, "c2", "D", "C", orientation={"roll": 20})
    member(d, "b1", "B", "M", sec="B")
    member(d, "b2", "M", "C", sec="B")
    member(d, "arm", "S1", "T")
    member(d, "strut", "S2", "D", type="truss")
    d["constraints"] = {"link": {"type": "rigidLink", "master": "M", "slaves": ["S1", "S2"]}}
    case(d, "E", {"type": "nodal", "node": "S1", "Fx": 8, "Fy": -3, "Fz": -20, "Mx": 2, "My": -1, "Mz": 4},
         {"type": "nodal", "node": "T", "Fz": -5, "Fy": 6},
         {"type": "nodal", "node": "S2", "Fy": 7},
         {"type": "memberDistributed", "member": "arm", "direction": "Z", "w1": -3},
         {"type": "memberDistributed", "member": "b1", "direction": "Z", "w1": -8})
    return d, "E"


@model
def L08_diaphragm_z_up_unequal_columns():
    d = start()
    nodes(d, a=(0, 0, 0), b=(6, 0, -1), c=(6, 5, 0.5), e=(0, 5, 0), A=(0, 0, 4), B=(6, 0, 4), C=(6, 5, 4),
          E=(0, 5, 4), M=(2.2, 1.9, 4))
    for n in "abce":
        d["supports"][n] = {"type": "fixed" if n in "ab" else "pinned"}
    for k, (lo, up) in enumerate(zip("abce", "ABCE")):
        member(d, f"col{k}", lo, up, sec="B" if k == 2 else "S", orientation={"roll": 90 if k % 2 else 0})
    for k, (i, j) in enumerate((("A", "B"), ("B", "C"), ("C", "E"), ("E", "A"))):
        member(d, f"bm{k}", i, j)
    d["constraints"] = {"floor": {"type": "rigidDiaphragm", "master": "M", "slaves": ["A", "B", "C", "E"]}}
    d["supports"]["M"] = {"dofs": {"uz": "fixed", "rx": "fixed", "ry": "fixed"}}
    case(d, "H", {"type": "nodal", "node": "M", "Fx": 50, "Fy": 20, "Mz": 30},
         {"type": "nodal", "node": "C", "Fx": 15, "Fz": -40},
         {"type": "memberDistributed", "member": "bm0", "direction": "Z", "w1": -10},
         {"type": "memberDistributed", "member": "bm1", "direction": "X", "w1": 4})
    return d, "H"


@model
def L09_two_storey_diaphragms_y_up():
    d = start(("kN", "mm"), "Y")
    xz = [(0, 0), (7000, 500), (6500, 6000), (-500, 5000)]
    for k, (x, z) in enumerate(xz):
        nodes(d, **{f"g{k}": (x, 0, z), f"p{k}": (x, 3500, z), f"q{k}": (x, 6800, z)})
        d["supports"][f"g{k}"] = {"type": "fixed"}
        member(d, f"ca{k}", f"g{k}", f"p{k}", orientation={"roll": 45 * k})
        member(d, f"cb{k}", f"p{k}", f"q{k}")
    for lev in "pq":
        for k in range(4):
            member(d, f"b{lev}{k}", f"{lev}{k}", f"{lev}{(k + 1) % 4}", sec="B")
    d["constraints"] = {"f1": {"type": "rigidDiaphragm", "master": "p0", "slaves": ["p1", "p2", "p3"]},
                        "f2": {"type": "rigidDiaphragm", "master": "q2", "slaves": ["q0", "q1", "q3"], "normal": "Y"}}
    case(d, "E", {"type": "nodal", "node": "p1", "Fx": 80, "Fz": -30}, {"type": "nodal", "node": "q3", "Fx": 120, "Fz": 45},
         {"type": "nodal", "node": "q0", "My": 5.0e4},
         *[{"type": "memberDistributed", "member": f"bq{k}", "direction": "Y", "w1": -0.02} for k in range(4)])
    return d, "E"


@model
def L10_diaphragm_normal_x_wall():
    d = start()
    nodes(d, a=(0, 0, 0), b=(0, 4, 0), A=(0, 0, 3), B=(0, 4, 3), C=(0, 2, 5), P=(5, 0, 3), Q=(5, 4, 3), p=(5, 0, 0), q=(5, 4, 0))
    for n in "abpq":
        d["supports"][n] = {"type": "fixed"}
    member(d, "w1", "a", "A")
    member(d, "w2", "b", "B")
    member(d, "w3", "A", "C")
    member(d, "w4", "B", "C")
    member(d, "c1", "p", "P")
    member(d, "c2", "q", "Q")
    member(d, "t1", "A", "P", sec="B")
    member(d, "t2", "B", "Q", sec="B")
    d["constraints"] = {"wall": {"type": "rigidDiaphragm", "master": "A", "slaves": ["B", "C"], "normal": "X"}}
    case(d, "E", {"type": "nodal", "node": "C", "Fy": 30, "Fz": -10, "Fx": 4, "Mx": 12},
         {"type": "nodal", "node": "B", "Fy": -8},
         {"type": "memberDistributed", "member": "t1", "direction": "Z", "w1": -6})
    return d, "E"


@model
def L11_support_springs_stiff_and_soft():
    d = start(("N", "mm"))
    nodes(d, A=(0, 0, 0), B=(4000, 0, 0), C=(8000, 1000, 500), D=(12000, 0, 0), T=(8000, 1000, 3500))
    d["supports"] = {"A": {"dofs": {"ux": "fixed", "uy": "fixed", "uz": {"spring": 1.0e15}, "rx": "fixed",
                                    "ry": {"spring": 1.0e-9}, "rz": {"spring": 3.0e9}}},
                     "B": {"dofs": {"uz": {"spring": 2500.0}, "uy": {"spring": 1.0e-6}}},
                     "D": {"type": "pinned", "dofs": {"rx": {"spring": 1.0e18}}},
                     "T": {"dofs": {"ux": {"spring": 800.0}, "uy": {"spring": 1.0e12}, "rz": {"spring": 5.0e8}}}}
    member(d, "s1", "A", "B")
    member(d, "s2", "B", "C", orientation={"roll": 35})
    member(d, "s3", "C", "D")
    member(d, "post", "C", "T")
    case(d, "Q", {"type": "memberDistributed", "member": "s1", "direction": "Z", "w1": -12},
         {"type": "memberDistributed", "member": "s2", "direction": "Z", "w1": -12, "w2": -4},
         {"type": "nodal", "node": "T", "Fx": 15000, "Fy": 4000, "Mz": 2.0e6},
         {"type": "nodal", "node": "C", "Fy": 9000})
    return d, "Q"


@model
def L12_prescribed_support_displacements():
    d = start(("kip", "in"), "Y")
    nodes(d, A=(0, 0, 0), B=(0, 150, 0), C=(200, 170, 40), D=(210, 0, 50), F=(400, 150, 0), H=(400, 0, 0))
    d["supports"] = {"A": {"type": "fixed"}, "D": {"type": "pinned"}, "H": {"type": "fixed"}}
    member(d, "c1", "A", "B")
    member(d, "c2", "D", "C", orientation={"roll": 30})
    member(d, "c3", "H", "F")
    member(d, "b1", "B", "C", sec="B")
    member(d, "b2", "C", "F", sec="B", releases={"j": {"rz": True}})
    case(d, "S", {"type": "supportDisplacement", "node": "A", "uy": -0.75, "rz": 0.004, "ux": 0.2},
         {"type": "supportDisplacement", "node": "D", "uy": -1.1, "uz": 0.3},
         {"type": "supportDisplacement", "node": "H", "rx": -0.002, "ry": 0.003},
         {"type": "memberDistributed", "member": "b1", "direction": "Y", "w1": -0.1})
    return d, "S"


@model
def L13_trapezoid_partial_projected_on_inclined():
    d = start()
    nodes(d, A=(0, 0, 0), B=(0, 0, 3), R=(5, 1, 6), C=(10, 0, 3.5), D=(10, 0, 0), P=(5, 7, 5))
    d["supports"] = {"A": {"type": "fixed"}, "D": {"type": "fixed"}, "P": {"type": "pinned"}}
    member(d, "c1", "A", "B")
    member(d, "c2", "D", "C")
    member(d, "r1", "B", "R", orientation={"roll": 12})
    member(d, "r2", "R", "C")
    member(d, "pur", "R", "P", orientation={"referenceVector": [1, 0, 0.3]})
    case(d, "SN",
         {"type": "memberDistributed", "member": "r1", "direction": "Z", "w1": -6, "w2": -2, "projected": True},
         {"type": "memberDistributed", "member": "r2", "direction": "Z", "w1": -1, "w2": -7, "start": 0.15, "end": 0.85,
          "projected": True},
         {"type": "memberDistributed", "member": "r2", "direction": "X", "w1": 3, "w2": 0, "start": 1.0, "end": 4.0,
          "positionType": "absolute", "projected": True},
         {"type": "memberDistributed", "member": "pur", "direction": "Y", "w1": 2, "w2": 5, "end": 0.5},
         {"type": "memberDistributed", "member": "pur", "direction": "z", "w1": 1.5, "start": 0.5},
         {"type": "memberDistributed", "member": "c1", "direction": "x", "w1": -4, "w2": -1},
         {"type": "memberPoint", "member": "r1", "direction": "X", "value": 5, "position": 2.0, "positionType": "absolute"},
         {"type": "memberPoint", "member": "pur", "direction": "Z", "kind": "moment", "value": 3, "position": 0.7},
         {"type": "memberPoint", "member": "r2", "direction": "y", "kind": "moment", "value": -2, "position": 0.25})
    return d, "SN"


@model
def L14_self_weight_direction_subset_lbf_ft():
    d = start(("lbf", "ft"), "Y")
    nodes(d, A=(0, 0, 0), B=(3, 14, 2), C=(22, 18, 5), D=(25, 0, 0), F=(12, 16, 20), H=(12, 0, 24))
    for n in "ADH":
        d["supports"][n] = {"type": "fixed"}
    member(d, "c1", "A", "B")
    member(d, "c2", "D", "C", orientation={"roll": 50})
    member(d, "c3", "H", "F")
    member(d, "b1", "B", "C", sec="B")
    member(d, "b2", "C", "F", sec="B", offsets={"i": [-1.0, 0, 1.0]})
    member(d, "t", "F", "B", type="truss")
    d["loadCases"]["D"] = {"loads": [{"type": "selfWeight"}]}
    d["loadCases"]["Q"] = {"loads": [{"type": "selfWeight", "factor": 0.3, "direction": [1, 0, 0.5],
                                      "members": ["b1", "b2", "t"]}]}
    d["combinations"]["U"] = {"factors": {"D": 1.2, "Q": 1.6}}
    return d, "U"


@model
def L15_shear_deformation_deep_short_members():
    d = start(shear=True)
    d["sections"]["V"] = {"type": "properties", "A": 0.3, "Iy": 0.004, "Iz": 0.016, "J": 0.01, "Asy": 0.25, "Asz": 0.2}
    nodes(d, A=(0, 0, 0), B=(0, 0, 2), C=(1.5, 0, 2), D=(1.5, 0, 0), F=(1.5, 2, 2), H=(1.5, 2, 0))
    for n in "ADH":
        d["supports"][n] = {"type": "fixed"}
    member(d, "c1", "A", "B", sec="V")
    member(d, "c2", "D", "C", sec="V")
    member(d, "c3", "H", "F", sec="V", orientation={"roll": 90})
    member(d, "b1", "B", "C", sec="V")
    member(d, "b2", "C", "F", sec="V")
    case(d, "V", {"type": "memberDistributed", "member": "b1", "direction": "Z", "w1": -300},
         {"type": "memberDistributed", "member": "b2", "direction": "Z", "w1": -200, "start": 0.25, "end": 0.75},
         {"type": "memberPoint", "member": "c1", "direction": "X", "value": 400, "position": 0.5},
         {"type": "nodal", "node": "F", "Fx": 500, "Fy": 300, "Mz": 100})
    return d, "V"


@model
def L16_space_truss_with_member_loads():
    d = start(("kN", "mm"))
    nodes(d, a=(0, 0, 0), b=(3000, 0, 0), c=(1500, 2600, 0), A=(200, 100, 4000), B=(2800, 100, 4000),
          C=(1500, 2400, 4000), T=(1500, 900, 7000))
    for n in "abc":
        d["supports"][n] = {"type": "pinned"}
    pairs = [("a", "A"), ("b", "B"), ("c", "C"), ("a", "B"), ("b", "C"), ("c", "A"), ("A", "B"), ("B", "C"),
             ("C", "A"), ("A", "T"), ("B", "T"), ("C", "T")]
    for k, (i, j) in enumerate(pairs):
        member(d, f"t{k}", i, j, type="truss")
    case(d, "W", {"type": "selfWeight", "factor": 50.0}, {"type": "nodal", "node": "T", "Fx": 40, "Fy": -25, "Fz": -90},
         {"type": "memberDistributed", "member": "t6", "direction": "Z", "w1": -0.01, "w2": -0.03},
         {"type": "memberPoint", "member": "t10", "direction": "Y", "value": 12, "position": 0.3})
    return d, "W"


@model
def L17_grillage_torsion_point_moments():
    d = start()
    nodes(d, n1=(0, 0, 0), n2=(4, 0, 0), n3=(8, 0, 0), n4=(0, 3, 0), n5=(4, 3, 0), n6=(8, 3, 0), n7=(4, 6, 0), n8=(4, -2, 0))
    for n in ("n1", "n3", "n4", "n6", "n7"):
        d["supports"][n] = {"type": "fixed" if n in ("n1", "n7") else "pinned"}
    for k, (i, j) in enumerate((("n1", "n2"), ("n2", "n3"), ("n4", "n5"), ("n5", "n6"), ("n2", "n5"), ("n5", "n7"),
                                ("n8", "n2"))):
        member(d, f"g{k}", i, j, orientation={"roll": 0 if k % 2 else 0})
    case(d, "T", {"type": "memberPoint", "member": "g4", "direction": "x", "kind": "moment", "value": 15, "position": 0.4},
         {"type": "memberPoint", "member": "g0", "direction": "Y", "kind": "moment", "value": -8, "position": 0.5},
         {"type": "memberPoint", "member": "g3", "direction": "Z", "value": -30, "position": 0.25},
         {"type": "nodal", "node": "n8", "Fz": -12, "Mx": 5, "My": 3},
         {"type": "memberDistributed", "member": "g5", "direction": "Z", "w1": -9, "w2": -3})
    return d, "T"


@model
def L18_gerber_beam_with_hinges():
    d = start(("kN", "mm"), "Y")
    xs = [0, 5000, 6200, 11000, 12500, 18000]
    for k, x in enumerate(xs):
        nodes(d, **{f"n{k}": (x, 0, 0)})
    d["supports"] = {"n0": {"type": "fixed"}, "n1": {"type": "roller"}, "n3": {"dofs": {"uy": "fixed", "uz": "fixed"}},
                     "n5": {"dofs": {"uy": "fixed", "uz": "fixed", "rx": "fixed"}}}
    member(d, "s0", "n0", "n1")
    member(d, "s1", "n1", "n2", releases={"j": {"rz": True}})
    member(d, "s2", "n2", "n3")
    member(d, "s3", "n3", "n4")
    member(d, "s4", "n4", "n5", releases={"i": {"rz": True, "ry": True}})
    case(d, "Q", *[{"type": "memberDistributed", "member": f"s{k}", "direction": "Y", "w1": -0.03 - 0.005 * k}
                   for k in range(5)],
         {"type": "memberPoint", "member": "s2", "direction": "Y", "value": -60, "position": 0.5},
         {"type": "memberPoint", "member": "s4", "direction": "Z", "value": 10, "position": 0.3},
         {"type": "nodal", "node": "n2", "Fx": 25})
    return d, "Q"


@model
def L19_very_short_and_very_long_members():
    d = start()
    nodes(d, A=(0, 0, 0), B=(0, 0, 0.02), C=(60, 0, 0.02), D=(60, 0.015, 0.02), F=(60, 0.015, -8), G=(30, 0, 0.02))
    d["supports"] = {"A": {"type": "fixed"}, "F": {"type": "fixed"}}
    member(d, "stub1", "A", "B")
    member(d, "long1", "B", "G")
    member(d, "long2", "G", "C")
    member(d, "stub2", "C", "D")
    member(d, "col", "D", "F")
    case(d, "Q", {"type": "memberDistributed", "member": "long1", "direction": "Z", "w1": -2},
         {"type": "memberDistributed", "member": "long2", "direction": "Z", "w1": -2},
         {"type": "nodal", "node": "G", "Fy": 3, "Mx": 1}, {"type": "memberPoint", "member": "stub1", "direction": "X",
                                                              "value": 50, "position": 0.5})
    return d, "Q"


@model
def L20_zero_loads():
    d, _ = L05_offsets_off_axis_skew()
    d["loadCases"] = {"none": {"loads": []}, "zero": {"loads": [{"type": "nodal", "node": "C", "Fx": 0.0},
                                                               {"type": "selfWeight", "factor": 0.0}]}}
    d["combinations"] = {"Z": {"factors": {"none": 1.0, "zero": 1.5}}}
    return d, "Z"


@model
def L21_offsets_and_loads_on_constrained_slaves():
    d = start(("kip", "in"))
    nodes(d, a=(0, 0, 0), b=(300, 0, 0), c=(300, 240, 0), e=(0, 240, 0), A=(0, 0, 150), B=(300, 0, 150),
          C=(300, 240, 160), E=(0, 240, 150), M=(120, 100, 150), X=(150, 120, 190))
    for n in "abce":
        d["supports"][n] = {"type": "fixed"}
    for k, (lo, up) in enumerate(zip("abce", "ABCE")):
        member(d, f"col{k}", lo, up, offsets={"j": [0, 0, -9.0]})
    member(d, "b0", "A", "B", sec="B", offsets={"i": [7, 0, 0], "j": [-7, 0, 2]}, releases={"j": {"rz": True}})
    member(d, "b1", "B", "C", sec="B", offsets={"axes": "local", "i": [7, 1, 0], "j": [-7, 0, 0]})
    member(d, "b2", "C", "E", sec="B", offsets={"i": [-7, 0, 0]})
    member(d, "b3", "E", "A", sec="B")
    member(d, "mast", "M", "X", offsets={"i": [0, 0, 6]})
    d["constraints"] = {"floor": {"type": "rigidDiaphragm", "master": "M", "slaves": ["A", "B", "E"]},
                        "lnk": {"type": "rigidLink", "master": "C", "slaves": ["X"]}}
    d["supports"]["C"] = {"dofs": {"ux": {"spring": 40.0}}}
    case(d, "E", {"type": "nodal", "node": "A", "Fx": 20, "Fy": 9, "Mz": 300, "Fz": -30},
         {"type": "nodal", "node": "E", "Fx": -6, "Mx": 100}, {"type": "nodal", "node": "X", "Fy": 14, "Fz": -8, "My": 60},
         {"type": "memberDistributed", "member": "b0", "direction": "Z", "w1": -0.2},
         {"type": "memberDistributed", "member": "b1", "direction": "Z", "w1": -0.15})
    return d, "E"


@model
def L22_mixed_n_mm_y_up_no_shear_area_on_one_axis():
    d = start(("N", "mm"), "Y", shear=True)
    d["sections"]["S"].pop("Asz")
    nodes(d, A=(0, 0, 0), B=(500, 3000, 200), C=(4500, 3400, 0), D=(5000, 0, -300), F=(2500, 5000, 1500))
    d["supports"] = {"A": {"type": "pinned", "dofs": {"ry": "fixed"}}, "D": {"type": "fixed"},
                     "F": {"dofs": {"uz": {"spring": 900.0}, "ux": "fixed"}}}
    member(d, "c1", "A", "B", orientation={"roll": -20})
    member(d, "c2", "D", "C", sec="B")
    member(d, "b", "B", "C")
    member(d, "r1", "B", "F", orientation={"referenceVector": [0, 0, 1]})
    member(d, "r2", "F", "C", releases={"i": {"rz": {"spring": 2.0e8}}})
    case(d, "U", {"type": "selfWeight"}, {"type": "memberDistributed", "member": "b", "direction": "Y", "w1": -8, "w2": -14},
         {"type": "memberPoint", "member": "r1", "direction": "z", "value": 3000, "position": 0.5},
         {"type": "nodal", "node": "F", "Fx": 1000, "Fy": -5000})
    return d, "U"


@model
def L23_prescribed_displacement_of_a_rigid_link_master():
    d = start()
    nodes(d, M=(0, 0, 0), S1=(1.5, 0.5, 0.4), S2=(-1.0, 1.2, 0.0), T1=(5, 1, 3), T2=(-4, 3, 2.5), G1=(5, 1, 0), G2=(-4, 3, 0))
    d["supports"] = {"M": {"type": "fixed"}, "G1": {"type": "pinned"}, "G2": {"type": "fixed"}}
    member(d, "m1", "S1", "T1")
    member(d, "m2", "S2", "T2", orientation={"roll": 15})
    member(d, "c1", "G1", "T1")
    member(d, "c2", "G2", "T2")
    member(d, "tie", "T1", "T2")
    d["constraints"] = {"cap": {"type": "rigidLink", "master": "M", "slaves": ["S1", "S2"]}}
    case(d, "S", {"type": "supportDisplacement", "node": "M", "uz": -0.02, "rx": 0.004, "rz": -0.003, "ux": 0.01},
         {"type": "nodal", "node": "S1", "Fy": 12}, {"type": "memberDistributed", "member": "tie", "direction": "Z", "w1": -5})
    return d, "S"


@model
def L24_building_diaphragm_braces_hinged_beams_kip_ft():
    d = start(("kip", "ft"))
    pts = [(0, 0), (24, 0), (24, 18), (0, 18)]
    for k, (x, y) in enumerate(pts):
        nodes(d, **{f"g{k}": (x, y, 0), f"u{k}": (x, y, 13)})
        d["supports"][f"g{k}"] = {"type": "pinned"}
        member(d, f"c{k}", f"g{k}", f"u{k}", orientation={"roll": 90 * (k % 2)})
    nodes(d, cm=(11, 10, 13))
    d["supports"]["cm"] = {"dofs": {"uz": "fixed", "rx": "fixed", "ry": "fixed"}}
    for k in range(4):
        member(d, f"b{k}", f"u{k}", f"u{(k + 1) % 4}", sec="B", releases={"i": {"rz": True, "ry": True},
                                                                       "j": {"rz": True, "ry": True}})
    member(d, "x0", "g0", "u1", type="truss")
    member(d, "x1", "g1", "u2", type="truss")
    member(d, "x2", "g3", "u2", type="truss")
    d["constraints"] = {"roof": {"type": "rigidDiaphragm", "master": "cm", "slaves": ["u0", "u1", "u2", "u3"]}}
    case(d, "W", {"type": "nodal", "node": "cm", "Fx": 30, "Fy": 12, "Mz": 80},
         *[{"type": "memberDistributed", "member": f"b{k}", "direction": "Z", "w1": -1.1} for k in range(4)])
    return d, "W"


@model
def L25_inclined_supports_emulated_by_releases_and_springs():
    d = start(("kN", "m"), "Y")
    nodes(d, A=(0, 0, 0), B=(2, 3, 1), C=(6, 4, -1), D=(9, 0, 0), K=(6, 0, -1), P=(2, 3, 5))
    d["supports"] = {"A": {"type": "fixed"}, "D": {"type": "roller", "dofs": {"uz": "fixed", "rx": {"spring": 100.0}}},
                     "K": {"type": "fixed"}, "P": {"type": "pinned", "dofs": {"rz": "fixed"}}}
    member(d, "a", "A", "B", orientation={"roll": 33})
    member(d, "b", "B", "C", releases={"i": {"ry": {"spring": 800.0}}, "j": {"rx": True}})
    member(d, "c", "C", "D", orientation={"referenceVector": [0, 0, 1]})
    member(d, "k", "K", "C", releases={"j": {"ux": {"spring": 5000.0}, "rz": True}})
    member(d, "p", "P", "B", releases={"j": {"rx": True, "uy": {"spring": 300.0}, "uz": {"spring": 200.0}}})
    case(d, "Q", {"type": "memberDistributed", "member": "b", "direction": "Y", "w1": -9, "projected": True},
         {"type": "memberDistributed", "member": "c", "direction": "Z", "w1": 2, "w2": 6},
         {"type": "memberPoint", "member": "p", "direction": "Y", "value": -11, "position": 0.4},
         {"type": "nodal", "node": "C", "Fx": 7, "Fz": -4, "My": 2})
    return d, "Q"


@model
def L26_shear_deformation_uniform_loads_only():
    d = start(shear=True)
    d["sections"]["V"] = {"type": "properties", "A": 0.3, "Iy": 0.004, "Iz": 0.016, "J": 0.01, "Asy": 0.25, "Asz": 0.2}
    nodes(d, A=(0, 0, 0), B=(0.4, 0.2, 2), C=(2.0, 0.5, 2.3), D=(2.2, 0, 0), F=(1.8, 2.5, 2), H=(1.5, 2.6, 0))
    d["supports"] = {"A": {"type": "fixed"}, "D": {"type": "pinned"}, "H": {"type": "fixed"}}
    member(d, "c1", "A", "B", sec="V", orientation={"roll": 20})
    member(d, "c2", "D", "C", sec="V")
    member(d, "c3", "H", "F", sec="V", orientation={"roll": 90})
    member(d, "b1", "B", "C", sec="V", orientation={"roll": -15})
    member(d, "b2", "C", "F", sec="V")
    case(d, "V", {"type": "memberDistributed", "member": "b1", "direction": "Z", "w1": -300},
         {"type": "memberDistributed", "member": "b2", "direction": "y", "w1": -200},
         {"type": "memberDistributed", "member": "c1", "direction": "X", "w1": 80},
         {"type": "nodal", "node": "F", "Fx": 500, "Fy": 300, "Mz": 100}, {"type": "selfWeight", "factor": 40})
    return d, "V"


@model
def L27_skew_frame_bending_releases_and_offsets_uniform_loads():
    d = start(("kip", "in"))
    nodes(d, A=(0, 0, 0), B=(20, 10, 150), C=(260, 40, 170), D=(280, 0, 0), F=(150, 250, 160), H=(140, 260, 0))
    for n in "ADH":
        d["supports"][n] = {"type": "fixed"}
    member(d, "c1", "A", "B", orientation={"roll": 10})
    member(d, "c2", "D", "C")
    member(d, "c3", "H", "F", orientation={"roll": 60})
    member(d, "b1", "B", "C", sec="B", offsets={"i": [6, 1, 0], "j": [-8, -1, 0.5]}, releases={"i": {"rz": True}})
    member(d, "b2", "C", "F", sec="B", orientation={"roll": 25}, releases={"i": {"ry": True}, "j": {"rz": True, "ry": True}})
    member(d, "b3", "F", "B", sec="B", offsets={"axes": "local", "i": [5, 0, 0]}, releases={"j": {"rz": True}})
    case(d, "U", {"type": "memberDistributed", "member": "b1", "direction": "Z", "w1": -0.2},
         {"type": "memberDistributed", "member": "b2", "direction": "Z", "w1": -0.15},
         {"type": "memberDistributed", "member": "b3", "direction": "y", "w1": -0.1},
         {"type": "memberDistributed", "member": "b3", "direction": "X", "w1": 0.03},
         {"type": "nodal", "node": "C", "Fx": 10, "Fy": -6, "Mz": 100})
    return d, "U"


# ====================================================================== second order
SECOND = {}


def second(f):
    SECOND[f.__name__] = f
    return f


def _so(d, **kw):
    d.setdefault("analysis", {}).update({"convergenceTolerance": 1.0e-11, "maxIterations": 200, **kw})
    for c in d["combinations"].values():
        c["analysisType"] = "secondOrder"
    return d


@second
def S01_cantilever_beam_column_both_planes():
    d = start()
    nodes(d, A=(0, 0, 0), B=(0, 0, 4))
    d["supports"] = {"A": {"type": "fixed"}}
    member(d, "col", "A", "B", orientation={"roll": 30})
    case(d, "P", {"type": "nodal", "node": "B", "Fz": -450, "Fx": 12, "Fy": -7, "Mx": 3},
         {"type": "memberDistributed", "member": "col", "direction": "X", "w1": 4, "w2": 1},
         {"type": "memberPoint", "member": "col", "direction": "Y", "value": 6, "position": 0.6})
    return _so(d), "P"


@second
def S02_cantilever_beam_column_with_shear_deformation_and_tension():
    d = start(shear=True)
    d["sections"]["S"].update({"Asy": 0.0006, "Asz": 0.0009})       # weak in shear: N/GAs is not negligible
    nodes(d, A=(0, 0, 0), B=(0, 0, 4), C=(3, 0, 0), D=(3, 0, 4))
    d["supports"] = {"A": {"type": "fixed"}, "C": {"type": "fixed"}}
    member(d, "col", "A", "B")
    member(d, "tie", "C", "D")
    case(d, "P", {"type": "nodal", "node": "B", "Fz": -500, "Fx": 12, "Fy": -7},
         {"type": "nodal", "node": "D", "Fz": 2500, "Fx": 30, "Fy": 20},
         {"type": "memberDistributed", "member": "tie", "direction": "X", "w1": 9},
         {"type": "memberDistributed", "member": "col", "direction": "Y", "w1": 2, "start": 0.2, "end": 0.7})
    return _so(d), "P"


@second
def S03_portal_with_leaning_column():
    d = start()
    nodes(d, a=(0, 0, 0), A=(0, 0, 5), b=(8, 0, 0), B=(8, 0, 5), c=(14, 0, 0), C=(14, 0, 5))
    d["supports"] = {"a": {"type": "fixed"}, "b": {"type": "pinned", "dofs": {"rx": "fixed", "rz": "fixed"}},
                     "c": {"type": "pinned", "dofs": {"rz": "fixed"}}}
    member(d, "c1", "a", "A")
    member(d, "c2", "b", "B")
    member(d, "lean", "c", "C", releases={"j": {"ry": True, "rz": True}})
    member(d, "g1", "A", "B", sec="B")
    member(d, "g2", "B", "C", sec="B")
    for n in "ABC":
        d["supports"][n] = {"dofs": {"uy": "fixed"}}
    case(d, "U", {"type": "nodal", "node": "A", "Fx": 25, "Fz": -300}, {"type": "nodal", "node": "B", "Fz": -400},
         {"type": "nodal", "node": "C", "Fz": -900},
         {"type": "memberDistributed", "member": "g1", "direction": "Z", "w1": -20})
    return _so(d), "U"


@second
def S04_skew_space_frame_y_up_kip_in():
    d, cid = L02_space_frame_y_up_kip_in()
    d["loadCases"]["G"]["loads"] += [{"type": "nodal", "node": n, "Fy": -250.0} for n in "ABCE"]
    return _so(d), cid


@second
def S05_offsets_along_the_axis_carrying_compression():
    d = start()
    nodes(d, a=(0, 0, 0), A=(0, 0, 5), b=(6, 0, 0), B=(6, 0, 5))
    d["supports"] = {"a": {"type": "fixed"}, "b": {"type": "fixed"}}
    member(d, "c1", "a", "A", offsets={"i": [0, 0, 0.6], "j": [0, 0, -0.4]})
    member(d, "c2", "b", "B", offsets={"j": [0, 0, -0.5]}, orientation={"roll": 90})
    member(d, "g", "A", "B", sec="B", offsets={"i": [0.3, 0, 0], "j": [-0.3, 0, 0]})
    case(d, "U", {"type": "nodal", "node": "A", "Fx": 30, "Fy": 8, "Fz": -700}, {"type": "nodal", "node": "B", "Fz": -500, "Fy": 5},
         {"type": "memberDistributed", "member": "g", "direction": "Z", "w1": -15})
    return _so(d), "U"


@second
def S06_truss_members_with_transverse_load_and_partial_releases():
    d = start()
    nodes(d, a=(0, 0, 0), A=(0, 0, 4), b=(5, 0, 0), B=(5, 0, 4), e=(2.5, 4, 0))
    d["supports"] = {"a": {"type": "fixed"}, "b": {"type": "fixed"}, "e": {"type": "pinned"}}
    member(d, "c1", "a", "A")
    member(d, "c2", "b", "B", releases={"i": {"rz": {"spring": 4000.0}, "ry": {"spring": 9000.0}}})
    member(d, "g", "A", "B", sec="B", releases={"i": {"rz": {"spring": 20000.0}}})
    member(d, "strut", "e", "A", type="truss")
    member(d, "strut2", "e", "B", type="truss")
    case(d, "U", {"type": "nodal", "node": "A", "Fx": 20, "Fy": -60, "Fz": -500}, {"type": "nodal", "node": "B", "Fz": -600, "Fy": -70},
         {"type": "memberDistributed", "member": "strut", "direction": "X", "w1": 3},
         {"type": "memberPoint", "member": "strut2", "direction": "Z", "value": -5, "position": 0.4})
    return _so(d), "U"


@second
def S07_columns_under_self_weight_varying_axial_force():
    d = start(("kN", "mm"), "Y")
    nodes(d, a=(0, 0, 0), A=(0, 9000, 0), b=(5000, 0, 0), B=(5000, 9000, 0))
    d["supports"] = {"a": {"type": "fixed"}, "b": {"type": "fixed"}}
    member(d, "c1", "a", "A")
    member(d, "c2", "b", "B", orientation={"roll": 45})
    member(d, "g", "A", "B")
    case(d, "U", {"type": "selfWeight", "factor": 25.0}, {"type": "nodal", "node": "A", "Fx": 6, "Fz": 4, "Fy": -60},
         {"type": "memberDistributed", "member": "c2", "direction": "Y", "w1": -0.03, "w2": -0.01})
    return _so(d, convergenceTolerance=1.0e-9), "U"


@second
def S08_diaphragm_storey_with_leaning_columns():
    d = start()
    pts = [(0, 0), (7, 0), (7, 6), (0, 6), (3.5, 3)]
    for k, (x, y) in enumerate(pts):
        nodes(d, **{f"g{k}": (x, y, 0), f"u{k}": (x, y, 4.2)})
        d["supports"][f"g{k}"] = {"type": "fixed" if k < 2 else "pinned"}
        member(d, f"c{k}", f"g{k}", f"u{k}", orientation={"roll": 30 * k})
        if k >= 2:
            d["supports"][f"g{k}"] = {"type": "pinned", "dofs": {"rz": "fixed"}}
    d["constraints"] = {"floor": {"type": "rigidDiaphragm", "master": "u4", "slaves": ["u0", "u1", "u2", "u3"]}}
    case(d, "U", {"type": "nodal", "node": "u4", "Fx": 40, "Fy": 25, "Mz": 30, "Fz": -300},
         *[{"type": "nodal", "node": f"u{k}", "Fz": -120.0 - 20 * k} for k in range(4)])
    return _so(d), "U"


# ====================================================================== tension / compression only
ONE_SIDED = {}


def one_sided(f):
    ONE_SIDED[f.__name__] = f
    return f


@one_sided
def T01_two_storey_x_braced_frame():
    d = start()
    for k, z in enumerate((0, 3.5, 7.2)):
        nodes(d, **{f"l{k}": (0, 0, z), f"r{k}": (5, 0, z)})
    d["supports"] = {"l0": {"type": "fixed"}, "r0": {"type": "fixed"}}
    for k in (1, 2):
        member(d, f"cl{k}", f"l{k - 1}", f"l{k}")
        member(d, f"cr{k}", f"r{k - 1}", f"r{k}")
        member(d, f"g{k}", f"l{k}", f"r{k}", sec="B")
        member(d, f"x{k}a", f"l{k - 1}", f"r{k}", type="tensionOnly")
        member(d, f"x{k}b", f"r{k - 1}", f"l{k}", type="tensionOnly")
        for s in "lr":
            d["supports"][f"{s}{k}"] = {"dofs": {"uy": "fixed"}}
    case(d, "W", {"type": "nodal", "node": "l1", "Fx": 60}, {"type": "nodal", "node": "l2", "Fx": -35},
         {"type": "nodal", "node": "r2", "Fz": -200}, {"type": "memberDistributed", "member": "g1", "direction": "Z", "w1": -20})
    return d, "W"


@one_sided
def T02_space_tower_with_cross_bracing_y_up_kip_in():
    d = start(("kip", "in"), "Y")
    base = [(0, 0), (120, 0), (120, 100), (0, 100)]
    for lev, y in enumerate((0, 150, 290)):
        for k, (x, z) in enumerate(base):
            sx, sz = (1 - 0.08 * lev), (1 - 0.05 * lev)
            nodes(d, **{f"n{lev}{k}": (60 + (x - 60) * sx, y, 50 + (z - 50) * sz)})
    for k in range(4):
        d["supports"][f"n0{k}"] = {"type": "fixed"}
    for lev in (1, 2):
        for k in range(4):
            j = (k + 1) % 4
            member(d, f"leg{lev}{k}", f"n{lev - 1}{k}", f"n{lev}{k}")
            member(d, f"ring{lev}{k}", f"n{lev}{k}", f"n{lev}{j}")
            member(d, f"xa{lev}{k}", f"n{lev - 1}{k}", f"n{lev}{j}", type="tensionOnly")
            member(d, f"xb{lev}{k}", f"n{lev - 1}{j}", f"n{lev}{k}", type="tensionOnly")
    case(d, "W", {"type": "nodal", "node": "n20", "Fx": 25, "Fz": 12, "Fy": -40}, {"type": "nodal", "node": "n22", "Fx": 10, "Fz": -18},
         {"type": "nodal", "node": "n11", "Fz": 15, "My": 200})
    return d, "W"


@one_sided
def T03_compression_only_props_and_tension_only_ties_with_settlement():
    d = start(("kN", "mm"))
    nodes(d, A=(0, 0, 0), B=(4000, 0, 0), C=(8000, 0, 0), D=(12000, 0, 0), P1=(4000, 0, -3000), P2=(8000, 0, -3000),
          H1=(2000, 0, 3000), H2=(10000, 0, 3000))
    d["supports"] = {"A": {"type": "pinned", "dofs": {"rx": "fixed"}}, "D": {"type": "roller", "dofs": {"uy": "fixed"}},
                     "P1": {"type": "pinned"}, "P2": {"type": "pinned"}, "H1": {"type": "pinned"}, "H2": {"type": "pinned"}}
    member(d, "s1", "A", "B")
    member(d, "s2", "B", "C")
    member(d, "s3", "C", "D")
    member(d, "prop1", "P1", "B", type="compressionOnly")
    member(d, "prop2", "P2", "C", type="compressionOnly")
    member(d, "tie1", "H1", "B", type="tensionOnly")
    member(d, "tie2", "H2", "C", type="tensionOnly")
    case(d, "Q", {"type": "memberDistributed", "member": "s1", "direction": "Z", "w1": -0.03},
         {"type": "memberDistributed", "member": "s3", "direction": "Z", "w1": 0.012},
         {"type": "nodal", "node": "B", "Fx": 40}, {"type": "supportDisplacement", "node": "P1", "uz": -4.0})
    return d, "Q"
