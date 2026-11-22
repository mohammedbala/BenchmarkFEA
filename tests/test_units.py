"""Exact unit conversion; a file and its numerically converted twin give the same model."""
import copy
import re
from fractions import Fraction

import pytest

from benchfea import io
from benchfea.model import build_model
from benchfea.unitmap import convert_document, unit_map
from benchfea.units import FORCE_UNITS, INTERNAL, LENGTH_UNITS, UnitSystem, convert, parse_dimension

from helpers import assert_same

# Independent statement of the exact factors (NIST SP 811 App. B.8 and footnotes 22-23).
KIP_IN_N = Fraction("4448.2216152605")
IN_IN_MM = Fraction("25.4")
FT_IN_MM = Fraction("304.8")


def test_exact_base_factors():
    assert FORCE_UNITS["lbf"] == Fraction("4.4482216152605")
    assert FORCE_UNITS["lbf"] == Fraction("0.45359237") * Fraction("9.80665")
    assert FORCE_UNITS["kip"] == KIP_IN_N and FORCE_UNITS["kN"] == 1000 and FORCE_UNITS["N"] == 1
    assert LENGTH_UNITS["in"] == IN_IN_MM and LENGTH_UNITS["ft"] == FT_IN_MM == 12 * IN_IN_MM
    assert LENGTH_UNITS["m"] == 1000 and LENGTH_UNITS["mm"] == 1


def test_internal_system_is_newton_millimetre():
    assert (INTERNAL.force, INTERNAL.length) == ("N", "mm")
    for dim in ("force", "length", "force/length^2", "moment", "length^6"):
        assert INTERNAL.factor(dim) == 1


@pytest.mark.parametrize("dim,exp", [
    ("force", (1, 0)), ("length", (0, 1)), ("length^2", (0, 2)), ("length^4", (0, 4)), ("length^6", (0, 6)),
    ("force/length", (1, -1)), ("force/length^2", (1, -2)), ("force/length^3", (1, -3)), ("moment", (1, 1)),
    ("moment/rotation", (1, 1)), ("force*length", (1, 1)), ("dimensionless", (0, 0)), ("rotation", (0, 0)),
    ("degree", (0, 0)),
])
def test_parse_dimension(dim, exp):
    assert parse_dimension(dim) == exp


def test_parse_dimension_rejects_unknown():
    with pytest.raises(ValueError):
        parse_dimension("mass/length")


def test_derived_factors_are_exact_rationals():
    us = UnitSystem("kip", "in")
    assert us.factor("force/length^2") == KIP_IN_N / IN_IN_MM**2          # ksi -> MPa
    assert us.factor("moment") == KIP_IN_N * IN_IN_MM
    assert us.factor("length^4") == IN_IN_MM**4
    assert us.factor("force/length^3") == KIP_IN_N / IN_IN_MM**3
    si = UnitSystem("kN", "m")
    assert si.factor("force/length^2") == Fraction(1, 1000)                # kPa -> MPa
    assert si.factor("force/length") == 1                                  # kN/m == N/mm
    assert si.factor("moment") == 10**6
    assert UnitSystem("lbf", "ft").factor("moment") == Fraction("4.4482216152605") * FT_IN_MM


def test_known_values():
    us = UnitSystem("kip", "in")
    # 29000 ksi = 29000 * 4448.2216152605 / 645.16 MPa
    assert us.to_internal(29000.0, "force/length^2") == pytest.approx(199947.961501882478765, rel=1e-14)
    assert us.to_internal(1.0, "length") == 25.4
    assert us.from_internal(25.4, "length") == 1.0
    assert UnitSystem("kN", "mm").to_internal(210.0, "force/length^2") == 210000.0
    assert convert(1.0, "length", UnitSystem("kip", "ft"), UnitSystem("kip", "in")) == 12.0
    assert convert(3.0, "rotation", UnitSystem("kip", "ft"), UnitSystem("kN", "m")) == 3.0


@pytest.mark.parametrize("force", list(FORCE_UNITS))
@pytest.mark.parametrize("length", list(LENGTH_UNITS))
def test_round_trip_all_systems(force, length):
    us = UnitSystem(force, length)
    for dim in ("force", "length", "force/length^2", "moment", "length^4", "force/length^3", "moment/rotation"):
        for v in (1.0, 123.456, 2.9e4, 7.3e-5):
            assert us.from_internal(us.to_internal(v, dim), dim) == pytest.approx(v, rel=1e-15)


def test_unsupported_units_are_refused():
    with pytest.raises(ValueError):
        UnitSystem("kgf", "mm")
    with pytest.raises(ValueError):
        UnitSystem("kN", "cm")


def test_labels():
    us = UnitSystem("kip", "in")
    assert us.label("force/length^2") == "kip/in^2"
    assert us.label("moment") == "kip*in"
    assert us.label("moment/rotation") == "kip*in/rad"
    assert us.label("dimensionless") == "-"


# ----------------------------------------------------------------------------------------
# An independent, hand-written conversion of the fixture (no schema, no benchfea.units):
# every dimensional field of the format is listed here with its force/length exponents.
def hand_convert(doc, kf: float, kl: float, force: str, length: str):
    d = copy.deepcopy(doc)
    d["units"] = {"force": force, "length": length}

    def c(obj, key, f, l):
        if key in obj:
            obj[key] = obj[key] * kf**f * kl**l

    def spring(spec):
        for dof, v in spec.items():
            if isinstance(v, dict):
                c(v, "spring", 1, -1 if dof.startswith("u") else 1)

    for n in d["nodes"].values():
        for k in "xyz":
            c(n, k, 0, 1)
    for s in d.get("supports", {}).values():
        spring(s.get("dofs", {}))
    for m in d["materials"].values():
        for k in ("E", "G", "Fy", "Fu"):
            c(m, k, 1, -2)
        c(m, "unitWeight", 1, -3)
    powers = {"d": 1, "bf": 1, "tf": 1, "tw": 1, "r": 1, "H": 1, "B": 1, "t": 1, "D": 1, "h": 1, "b": 1,
              "A": 2, "Asy": 2, "Asz": 2, "Iy": 4, "Iz": 4, "J": 4, "Cw": 6, "Sy": 3, "Sz": 3, "Zy": 3, "Zz": 3}
    for s in d["sections"].values():
        for k, p in powers.items():
            c(s, k, 0, p)
            c(s.get("overrides", {}), k, 0, p)
    for m in d["members"].values():
        for end in ("i", "j"):
            spring(m.get("releases", {}).get(end, {}))
            if end in m.get("offsets", {}):
                m["offsets"][end] = [v * kl for v in m["offsets"][end]]
    for case in d["loadCases"].values():
        for ld in case["loads"]:
            absolute = ld.get("positionType") == "absolute"
            if ld["type"] == "nodal":
                for k in ("Fx", "Fy", "Fz"):
                    c(ld, k, 1, 0)
                for k in ("Mx", "My", "Mz"):
                    c(ld, k, 1, 1)
            elif ld["type"] == "memberPoint":
                c(ld, "value", 1, 1 if ld.get("kind") == "moment" else 0)
                if absolute:
                    c(ld, "position", 0, 1)
            elif ld["type"] == "memberDistributed":
                c(ld, "w1", 1, -1)
                c(ld, "w2", 1, -1)
                if absolute:
                    c(ld, "start", 0, 1)
                    c(ld, "end", 0, 1)
            elif ld["type"] == "supportDisplacement":
                for k in ("ux", "uy", "uz"):
                    c(ld, k, 0, 1)
    design = d.get("design", {})
    levels = [design.get("defaults", {})] + list(design.get("groups", {}).values()) \
        + list(design.get("members", {}).values())
    for lv in levels:
        for k in ("Lx", "Ly", "Lz", "Lb"):
            c(lv, k, 0, 1)
        c(lv, "An", 0, 2)
        for lim in lv.get("deflectionLimits", []):
            c(lim, "absolute", 0, 1)
        br = lv.get("bracing", {})
        if br.get("positionType") == "absolute":
            for k in ("major", "minor", "torsion", "topFlange", "bottomFlange"):
                if k in br:
                    br[k] = [v * kl for v in br[k]]
    return d


KIP_TO_KN = 4.4482216152605
IN_TO_MM = 25.4


def test_kip_in_file_and_converted_kn_mm_file_give_identical_models(fixtures):
    src = io.load(fixtures / "valid" / "frame_kip_in.json").data
    twin = hand_convert(src, KIP_TO_KN, IN_TO_MM, "kN", "mm")
    assert twin["nodes"]["N3"]["x"] == 240 * 25.4
    assert twin["materials"]["steel"]["E"] == pytest.approx(29000 * 4.4482216152605 / 645.16, rel=1e-14)
    a, b = build_model(src), build_model(twin)
    assert a.units == UnitSystem("kip", "in") and b.units == UnitSystem("kN", "mm")
    assert_same(a, b, ignore=("units", "defaults_used", "warnings", "source"))
    # the defaults-used reports agree once expressed in the same units
    assert [x.path for x in a.defaults_used] == [x.path for x in b.defaults_used]
    da = {x.path: x.value for x in a.defaults_used}
    db = {x.path: x.value for x in b.defaults_used}
    assert da["/design/members/C1/Lx"] == pytest.approx(144.0)
    assert db["/design/members/C1/Lx"] == pytest.approx(144.0 * 25.4)
    assert da["/design/members/C2/An"] == pytest.approx(10.0) and db["/design/members/C2/An"] == pytest.approx(6451.6)


def test_other_unit_systems_give_identical_models(fixtures):
    src = io.load(fixtures / "valid" / "frame_kip_in.json").data
    ref = build_model(src)
    for force, length, kf, kl in (("lbf", "ft", 1000.0, 1 / 12), ("N", "m", 4448.2216152605, 0.0254),
                                  ("kN", "m", 4.4482216152605, 0.0254), ("N", "mm", 4448.2216152605, 25.4),
                                  ("kip", "ft", 1.0, 1 / 12)):
        twin = hand_convert(src, kf, kl, force, length)
        assert_same(ref, build_model(twin), ignore=("units", "defaults_used", "warnings", "source"))


def test_schema_driven_conversion_agrees_with_hand_conversion(fixtures):
    src = io.load(fixtures / "valid" / "frame_kip_in.json").data
    auto = convert_document(src, UnitSystem("kN", "mm"))
    hand = hand_convert(src, KIP_TO_KN, IN_TO_MM, "kN", "mm")

    def walk(x, y, path=""):
        if isinstance(x, dict):
            assert x.keys() == y.keys(), path
            for k in x:
                walk(x[k], y[k], f"{path}/{k}")
        elif isinstance(x, list):
            assert len(x) == len(y), path
            for k, (p, q) in enumerate(zip(x, y)):
                walk(p, q, f"{path}/{k}")
        elif isinstance(x, (int, float)) and not isinstance(x, bool):
            assert x == pytest.approx(y, rel=1e-13, abs=0), path
        else:
            assert x == y, path

    walk(auto, hand)
    assert src["units"] == {"force": "kip", "length": "in"}        # the source is not modified


def test_unit_map_of_special_fields(fixtures):
    src = io.load(fixtures / "valid" / "frame_kip_in.json").data
    um = unit_map(src)
    assert um["/loadCases/L/loads/0/value"] == "force"
    assert um["/loadCases/L/loads/1/value"] == "moment"                 # kind = moment
    assert um["/loadCases/L/loads/0/position"] == "length"              # positionType absolute
    assert "/loadCases/L/loads/1/position" not in um                    # relative: a fraction
    assert "/loadCases/D/loads/3/start" not in um
    assert um["/loadCases/L/loads/2/end"] == "length"
    assert um["/supports/N4/dofs/ry/spring"] == "moment/rotation"
    assert um["/supports/N4/dofs/ux/spring"] == "force/length"
    assert um["/members/B1b/releases/j/rz/spring"] == "moment/rotation"
    assert um["/members/B1b/releases/j/ux/spring"] == "force/length"
    assert um["/design/groups/beams/bracing/minor/1"] == "length"       # absolute brace positions
    assert "/design/members/C1/bracing/torsion/0" not in um             # relative brace positions
    assert "/members/C2/orientation/roll" not in um                     # degrees
    assert "/loadCases/S/loads/0/ry" not in um                          # radians
    assert um["/loadCases/S/loads/0/uz"] == "length"
    assert um["/sections/beam/overrides/J"] == "length^4"
    assert not [p for p in um if re.search(r"/(factors|Kx|U|modes|resultStations)(/|$)", p)]


def test_numbers_are_converted_never_relabelled(fixtures):
    src = io.load(fixtures / "valid" / "beam.json").data
    relabelled = copy.deepcopy(src)
    relabelled["units"] = {"force": "N", "length": "mm"}                # same numbers, other label
    a, b = build_model(src), build_model(relabelled)
    assert a.materials["S355"].E == pytest.approx(2.1e5)                # 2.1e8 kN/m^2 = 2.1e5 MPa
    assert b.materials["S355"].E == pytest.approx(2.1e8)
    assert a.members["m1"].length == pytest.approx(6000.0) and b.members["m1"].length == pytest.approx(6.0)
    w = a.load_cases["G"].loads[0].w1
    assert abs(w).max() == pytest.approx(12.5)                          # kN/m == N/mm
