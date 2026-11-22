"""Solver-independent parts of benchfea.crosscheck: load preparation, fidelity records, the comparison
and the tier-C decision (with synthetic results), the capability table, the subprocess runner.
Only the runner tests at the end need a foreign package (they skip without it)."""
import importlib
import json
import math
import sys

import numpy as np
import pytest

from benchfea.crosscheck import (CAPABILITIES, APPROXIMATE, EXACT, UNSUPPORTED, Fidelity, ForeignResult, Tolerances,
                                 applied_resultant, available, capability, compare, cross_check, equilibrium,
                                 relative_difference, run_job)
from benchfea.crosscheck import capabilities as cap
from benchfea.crosscheck import opensees as ops_tr
from benchfea.crosscheck import pynite as pyn_tr
from benchfea.crosscheck import runner
from benchfea.crosscheck.prepare import prepare
from benchfea.model import build_model

from crosscheck_helpers import E, SEC, document, needs_opensees, needs_pynite


def _model(**kw):
    nodes = {"A": (0, 0, 0), "B": (0, 0, 3000), "C": (4000, 0, 3000), "D": (4000, 0, 0)}
    members = kw.pop("members", {"c1": ("A", "B"), "g": ("B", "C"), "c2": ("D", "C")})
    cases = kw.pop("cases", {
        "D": [{"type": "selfWeight"}, {"type": "memberDistributed", "member": "g", "direction": "Z", "w1": -4.0,
                                       "w2": -10.0, "start": 0.25, "end": 0.75}],
        "W": [{"type": "nodal", "node": "B", "Fx": 1e4}, {"type": "memberPoint", "member": "c1", "direction": "X",
                                                           "value": 2e3, "position": 0.5},
              {"type": "supportDisplacement", "node": "A", "uz": -1.0}]})
    return build_model(document(nodes, members, kw.pop("supports", {"A": "fixed", "D": "pinned"}), cases,
                                combos=kw.pop("combos", {"U": {"factors": {"D": 1.2, "W": 1.6}}}), **kw))


# ====================================================================================== prepare
def test_prepare_factors_and_local_loads():
    m = _model()
    p = prepare(m, "U")
    assert not p.second_order and p.alpha == 1.0
    assert p.load[1] == pytest.approx([1.6e4, 0, 0, 0, 0, 0])
    assert p.prescribed[0, 2] == pytest.approx(-1.6) and p.fixed[0].all()
    g = p.members["g"]                                       # beam along X: y = +Z
    w_self = 7.7e-5 * SEC["A"]
    dists = sorted(g.dists, key=lambda d: d[3] - d[2])
    (w1, w2, x1, x2), (s1, s2, a1, a2) = dists               # the partial load, then the self weight
    assert (x1, x2) == (1000.0, 3000.0) and w1[1] == pytest.approx(-4.8) and w2[1] == pytest.approx(-12.0)
    assert (a1, a2) == (0.0, 4000.0) and s1 == pytest.approx([0.0, -1.2 * w_self, 0.0])
    c1 = p.members["c1"]                                     # column: x = +Z, y = +X; gravity is axial
    assert c1.dists[0][0] == pytest.approx([-1.2 * w_self, 0.0, 0.0])
    kind, vec, x = c1.points[0]
    assert kind == "force" and vec == pytest.approx([0.0, 3.2e3, 0.0]) and x == 1500.0
    assert g.GAy is None and g.EIz == pytest.approx(E * SEC["Iz"])


def test_prepare_second_order_alpha_multiplies_all_loads():
    m = _model(combos={"U": {"factors": {"D": 1.0, "W": 1.0}, "analysisType": "secondOrder"}},
               analysis={"alpha": 1.6})
    p = prepare(m, "U")
    assert p.second_order and p.alpha == 1.6
    assert p.load[1, 0] == pytest.approx(1.6e4) and p.prescribed[0, 2] == pytest.approx(-1.6)
    R, _ = applied_resultant(m, "U")                         # resultant is reported for the un-multiplied loads
    assert R[0] == pytest.approx(1e4 + 2e3)


def test_prepare_stiffness_reduction_and_tau_b():
    kw = dict(combos={"U": {"factors": {"D": 1.0}, "analysisType": "secondOrder", "stiffnessReduction": True}},
              materials={"M": {"E": E, "G": 77000.0, "Fy": 345.0}})
    m = _model(analysis={"alpha": 1.0, "stiffnessReduction": {"factor": 0.8, "tauB": 0.9, "members": ["c1"]}}, **kw)
    p = prepare(m, "U")
    assert p.members["c1"].E == pytest.approx(0.8 * E) and p.members["c1"].EIz == pytest.approx(0.8 * 0.9 * E * SEC["Iz"])
    assert p.members["g"].E == pytest.approx(E) and p.members["g"].EIz == pytest.approx(E * SEC["Iz"])
    assert all(level == EXACT for f, level, _ in p.issues if f == "stiffnessReduction")
    # a "computed" tau_b depends on the analysis being checked: unsupported unless supplied
    m = _model(analysis={"alpha": 1.0, "stiffnessReduction": {"factor": 0.8, "tauB": "computed"}}, **kw)
    assert ("stiffnessReduction", UNSUPPORTED) in [(f, lv) for f, lv, _ in prepare(m, "U").issues]
    p = prepare(m, "U", tau_b={"c1": 0.7, "g": 1.0, "c2": 1.0})
    assert ("stiffnessReduction", UNSUPPORTED) not in [(f, lv) for f, lv, _ in p.issues]
    assert p.members["c1"].EIy == pytest.approx(0.8 * 0.7 * E * SEC["Iy"])


def test_prepare_notional_loads_and_one_way_members_are_not_guessed():
    m = _model(combos={"U": {"factors": {"D": 1.0}, "notionalLoad": {"direction": [1, 0, 0]}}})
    assert ("notionalLoad", UNSUPPORTED) in [(f, lv) for f, lv, _ in prepare(m, "U").issues]
    p = prepare(m, "U", extra_nodal_loads={"B": [50.0, 0, 0, 0, 0, 0]})
    assert ("notionalLoad", EXACT) in [(f, lv) for f, lv, _ in p.issues] and p.load[1, 0] == pytest.approx(50.0)
    m = _model(members={"c1": ("A", "B"), "g": ("B", "C"), "c2": ("D", "C"),
                        "x": {"i": "A", "j": "C", "type": "tensionOnly"}})
    assert ("memberTypes.tensionOnly", UNSUPPORTED) in [(f, lv) for f, lv, _ in prepare(m, "U").issues]
    for tr in (ops_tr, pyn_tr):
        assert tr.translate(m, "U").fidelity.level == UNSUPPORTED


def test_prepare_inactive_degrees_of_freedom():
    """Rotations of a node joined only by truss members, and the free DOFs of a memberless diaphragm
    master that drive nothing, are what benchfea reports as zero."""
    nodes = {"A": (0, 0, 0), "B": (4000, 0, 0), "C": (2000, 0, 1500), "T": (2000, 0, 4500), "M": (2000, 800, 4500)}
    members = {"t1": {"i": "A", "j": "C", "type": "truss"}, "t2": {"i": "B", "j": "C", "type": "truss"},
               "col": ("C", "T")}
    m = build_model(document(nodes, members, {"A": "pinned", "B": "pinned"},
                             {"L": [{"type": "nodal", "node": "M", "Fx": 1e3}]},
                             constraints={"d": {"type": "rigidDiaphragm", "master": "M", "slaves": ["T"]}}))
    p = prepare(m, "L")
    assert p.inactive == {"A": [3, 4, 5], "B": [3, 4, 5], "M": [2, 3, 4]}


def test_translation_jobs_are_plain_json_and_need_no_foreign_import():
    m = _model()
    before = set(sys.modules)
    for tr in (ops_tr, pyn_tr):
        t = tr.translate(m, "U")
        job = json.loads(json.dumps(t.job))
        assert job["solver"] == t.solver and isinstance(job["model"], list)
    assert not any(name.split(".")[0] in ("openseespy", "Pynite") for name in set(sys.modules) - before)
    cmds = ops_tr.translate(m, "U").job["model"]
    # OpenSees member axes: vecxz = benchfea local z
    transf = [c for c in cmds if c[0] == "geomTransf"]
    assert transf[0][1] == "Linear" and transf[0][3:6] == pytest.approx(list(m.members["c1"].R[2]))


def test_opensees_translation_choices():
    """Which OpenSees facilities a model is mapped to (no solver run)."""
    m = _model(members={"c1": ("A", "B"), "c2": ("D", "C"),
                        "g": {"i": "B", "j": "C", "releases": {"i": {"rz": True}, "j": {"rz": {"spring": 1e9}}}}})
    t = ops_tr.translate(m, "U")
    beams = [c for c in t.job["model"] if c[0] == "element" and c[1] == "elasticBeamColumn"]
    assert len(beams) == 2 + 2                               # g = core + one hinged partner
    assert sum("-releasez" in c for c in beams) == 1 and not t.controls
    loads = [c for c in t.job["pattern"] if c[0] == "eleLoad"]
    released = {c[2] for c in beams if c[3:5] == beams[2][3:5]}    # element tags between B and C
    assert all(c[2] not in released for c in loads)          # no eleLoad on a released member
    assert t.fidelity.level == EXACT
    t2 = ops_tr.translate(m, "U", native_loads=False)
    assert not [c for c in t2.job["pattern"] if c[0] == "eleLoad"]
    # second order: PDelta transformation, subdivision control, Newton
    m2 = _model(combos={"U": {"factors": {"D": 1.0}, "analysisType": "secondOrder"}}, analysis={"alpha": 1.0})
    t = ops_tr.translate(m2, "U")
    assert {c[1] for c in t.job["model"] if c[0] == "geomTransf"} == {"PDelta"}
    assert t.controls["subdivisions"] == {"coarse": 8, "fine": 16, "order": 2}
    assert len([c for c in t.job["model"] if c[0] == "element"]) == 3 * 16
    assert t.fidelity.items["secondOrder"].level == APPROXIMATE and t.fidelity.estimated_error is None


def test_pynite_translation_choices():
    m = _model(members={"c1": {"i": "A", "j": "B", "offsets": {"j": [0, 0, -200]}}, "c2": ("D", "C"),
                        "g": {"i": "B", "j": "C", "releases": {"i": {"rz": True, "ux": {"spring": 1e5}}}}},
               analysis={"shearDeformation": True})
    t = pyn_tr.translate(m, "U")
    feats = {i.feature: i.level for i in t.fidelity.items.values()}
    assert feats["shearDeformation"] == UNSUPPORTED and feats["releases.partial"] == UNSUPPORTED
    assert feats["offsets"] == APPROXIMATE and feats["releases.full"] == EXACT
    assert t.controls["stiffness"]["fine"] == 1e6 and t.controls["stiffness"]["order"] is None
    assert t.job["orient"] and all(len(y) == 3 for _, y in t.job["orient"])


def test_pynite_refuses_a_node_lying_on_a_member():
    """PyNite connects every node on a member's axis to that member; benchfea does not."""
    nodes = {"A": (0, 0, 0), "B": (4000, 0, 0), "X": (2000, 0, 0), "T": (2000, 0, 3000)}
    m = build_model(document(nodes, {"g": ("A", "B"), "c": ("X", "T")}, {"A": "fixed", "B": "fixed", "T": "fixed"},
                             {"L": [{"type": "nodal", "node": "X", "Fz": -1e3}]}))
    t = pyn_tr.translate(m, "L")
    assert "topology" in [i.feature for i in t.fidelity.unsupported]
    assert ops_tr.translate(m, "L").fidelity.level == EXACT


def test_fixture_models_translate_without_guessing(fixtures):
    """The repository fixtures use every feature of the format (one-way members, notional loads,
    stiffness reduction, second order, diaphragm with slaves off the master's plane, ...): both
    translators must produce a job for every combination and flag what they cannot represent."""
    for name in ("beam.json", "frame_kip_in.json"):
        model = build_model(fixtures / "valid" / name)
        for cid, combo in model.combinations.items():
            for tr in (ops_tr, pyn_tr):
                t = tr.translate(model, cid)
                json.dumps(t.job)
                feats = {i.feature: i.level for i in t.fidelity.items.values()}
                if any(mm.type in ("tensionOnly", "compressionOnly") for mm in model.members.values()):
                    assert feats["memberTypes.tensionOnly"] == UNSUPPORTED
                if combo.notional is not None:
                    assert feats["notionalLoad"] == UNSUPPORTED
                if combo.analysis_type == "secondOrder" and any(mm.has_offsets for mm in model.members.values()):
                    assert feats["offsets.secondOrder"] == APPROXIMATE and t.fidelity.estimated_error is None


# ====================================================================================== fidelity
def test_fidelity_record():
    f = Fidelity("x")
    f.add("a", EXACT, "one")
    f.add("b", APPROXIMATE, "two", control="penalty")
    assert f.level == APPROXIMATE and f.estimated_error is None
    f.set_error("penalty", 1e-5)
    f.add("c", APPROXIMATE, "three", error=2e-5)
    assert f.estimated_error == pytest.approx(3e-5)
    f.add("a", UNSUPPORTED, "worse")
    assert f.level == UNSUPPORTED and [i.feature for i in f.unsupported] == ["a"]
    assert "one" in f.items["a"].note and "worse" in f.items["a"].note
    g = Fidelity.from_dict(json.loads(json.dumps(f.to_dict())))
    assert g.to_dict() == f.to_dict()


# ====================================================================================== compare / tier C
def _result(solver, model, scale=1.0, level=EXACT, error=None, ok=True, shift=None):
    n = len(model.nodes)
    rng = np.random.default_rng(3)
    u = rng.normal(size=(n, 6)) * np.r_[[10.0] * 3, [0.01] * 3]
    r = np.zeros((n, 6))
    r[0] = rng.normal(size=6) * np.r_[[1e4] * 3, [1e7] * 3]
    ef = {mid: rng.normal(size=12) * np.r_[[1e4] * 3, [1e7] * 3, [1e4] * 3, [1e7] * 3] for mid in model.members}
    fid = Fidelity(solver)
    fid.add("members.frame", level, "synthetic", error=error)
    res = ForeignResult(solver, "U", ok, fid, u * scale, r * scale, {k: v * scale for k, v in ef.items()},
                        "" if ok else "boom")
    if shift:
        shift(res)
    return res


def test_compare_reports_largest_differences_per_quantity():
    m = _model()
    a = _result("a", m)

    def shift(res):
        res.displacements[2, 1] += 3e-4 * np.abs(res.displacements[:, :3]).max()
        res.end_forces["g"][11] *= 1.0 + 1e-3

    b = _result("b", m, shift=shift)
    c = compare(a, b, m)
    assert not c.passed and set(c.failures) == {"displacement", "endMoment"}
    assert c.worst["displacement"].where == "node C uy" and c.worst["displacement"].relative == pytest.approx(3e-4)
    assert c.worst["endMoment"].where == "member g end j Mz"
    assert c.largest[0].relative == c.max_relative and len(c.largest) == 2
    assert c.worst["rotation"].relative == 0.0 and "above tolerance" in c.summary()
    assert compare(a, b, m, Tolerances().widened(1e-2)).passed
    assert relative_difference(a, b, m) == c.max_relative
    with pytest.raises(ValueError):
        compare(a, _result("c", m, ok=False), m)


def test_compare_scales_do_not_amplify_noise():
    """A truss-like result with round-off moments: the moment scale has a floor tied to the force scale
    and the model size, so noise is not reported as a difference."""
    m = _model()
    a = _result("a", m)
    for res in (a,):
        res.reactions[:, 3:] = 0.0
        for v in res.end_forces.values():
            v[[3, 4, 5, 9, 10, 11]] = 0.0
    b = _result("b", m)
    b.reactions[:, 3:] = 1e-9
    for v in b.end_forces.values():
        v[[3, 4, 5, 9, 10, 11]] = 1e-7
    c = compare(a, b, m)
    assert c.passed and c.worst["endMoment"].relative < 1e-12


def test_tier_c_decision_rules():
    m = _model()
    exact = lambda name: _result(name, m)
    # two solvers represent the model and agree
    rep = cross_check(m, "U", solvers=("opensees", "pynite"), foreign={"opensees": exact("opensees"),
                                                                        "pynite": exact("pynite")})
    assert rep.cross_checked and rep.eligible == ["opensees", "pynite"]
    # only one solver: never cross-checked
    rep = cross_check(m, "U", solvers=("opensees",), foreign={"opensees": exact("opensees")})
    assert not rep.cross_checked and any("two are required" in r for r in rep.reasons)
    # an unsupported feature removes the solver even if its numbers agree
    rep = cross_check(m, "U", foreign={"opensees": exact("opensees"), "pynite": _result("pynite", m, level=UNSUPPORTED)})
    assert not rep.cross_checked and rep.eligible == ["opensees"]
    # an approximation without estimate does not count
    rep = cross_check(m, "U", foreign={"opensees": exact("opensees"), "pynite": _result("pynite", m, level=APPROXIMATE)})
    assert not rep.cross_checked and any("without an error estimate" in r for r in rep.reasons)
    # a converged approximation counts, with the tolerance widened by its estimate
    close = _result("pynite", m, scale=1.0 + 4e-4, level=APPROXIMATE, error=5e-4)
    rep = cross_check(m, "U", foreign={"opensees": exact("opensees"), "pynite": close})
    assert rep.cross_checked and rep.comparisons[0].tolerances.displacement == pytest.approx(1e-6 + 5e-4)
    # ... but not when the difference exceeds the estimate
    far = _result("pynite", m, scale=1.0 + 2e-3, level=APPROXIMATE, error=5e-4)
    rep = cross_check(m, "U", foreign={"opensees": exact("opensees"), "pynite": far})
    assert not rep.cross_checked and any("differ" in r for r in rep.reasons)
    # ... nor when the estimate itself is too large
    rep = cross_check(m, "U", foreign={"opensees": exact("opensees"),
                                       "pynite": _result("pynite", m, level=APPROXIMATE, error=2e-2)})
    assert not rep.cross_checked and any("exceeds" in r for r in rep.reasons)
    # a failed solver
    rep = cross_check(m, "U", foreign={"opensees": exact("opensees"), "pynite": _result("pynite", m, ok=False)})
    assert not rep.cross_checked and any("failed: boom" in r for r in rep.reasons)


def test_hand_equilibrium_from_load_definitions():
    """applied_resultant integrates the loads from their definition.  Gravity loads: Fz = -(self weight
    + trapezoid), Fx = wind; a result whose reactions balance them has zero residual."""
    m = _model()
    R, S = applied_resultant(m, "U")
    w = 7.7e-5 * SEC["A"]
    assert R[2] == pytest.approx(-1.2 * (w * (3000 + 4000 + 3000) + 0.5 * (4.0 + 10.0) * 2000.0))
    assert R[0] == pytest.approx(1.6 * (1e4 + 2e3))
    # moment about Y of the wind loads (z x Fx): node load at z = 3000, point load at z = 1500;
    # of the gravity loads (-x Fz): columns at x = 0 and 4000, beam loads at their centroids
    trap_x = 1000.0 + 2000.0 * (4.0 + 2 * 10.0) / (3 * (4.0 + 10.0))
    My = 1.6 * (1e4 * 3000 + 2e3 * 1500) + 1.2 * (w * 3000 * 4000 + w * 4000 * 2000 + 0.5 * 14.0 * 2000.0 * trap_x)
    assert R[4] == pytest.approx(My)
    res = _result("x", m)
    res.reactions[:] = 0.0
    res.reactions[0, :3] = -R[:3]
    res.reactions[0, 3:] = -R[3:]                            # node A is the origin
    eq = equilibrium(m, res)
    assert eq["force"] < 1e-14 and eq["moment"] < 1e-14
    res.reactions[0, 2] *= 1.01
    assert equilibrium(m, res)["force"] > 1e-3


# ====================================================================================== capabilities
def test_capability_table_is_complete_and_consistent():
    assert set(CAPABILITIES) == {"opensees", "pynite", "frame3dd"}
    for solver, table in CAPABILITIES.items():
        assert set(table) == set(cap.FEATURES), solver
        for feature, entry in table.items():
            assert entry["support"] in ("exact", "approximate", "unsupported", "not-run")
            assert entry["evidence"] in ("measured", "documentation", "source") and entry["how"]
    assert capability("pynite", "shearDeformation")["support"] == "unsupported"
    assert capability("opensees", "secondOrderWithShearDeformation")["support"] == "unsupported"
    assert capability("opensees")["rigidLink"]["support"] == "exact"
    assert all(v["support"] == "not-run" for v in CAPABILITIES["frame3dd"].values())
    assert cap.FORMULATION["frame3dd"]["runnable"] is False
    for urls in cap.SOURCES.values():
        assert urls and all(u.startswith("https://") for u in urls)
    json.dumps({"c": CAPABILITIES, "f": cap.FORMULATION, "s": cap.SOURCES})


# ====================================================================================== runner
def test_missing_package_is_reported_not_raised(monkeypatch):
    monkeypatch.setitem(runner.SOLVERS, "opensees", ("benchfea.crosscheck.opensees", "no_such_package_xyz"))
    assert not available("opensees")
    out = run_job({"solver": "opensees"})
    assert out["ok"] is False and out["missing"] and "not installed" in out["error"]
    m = _model()
    rep = cross_check(m, "U", solvers=("opensees",))
    assert not rep.cross_checked and "not installed" in rep.reasons[0]
    with pytest.raises(ValueError):
        runner.run("nastran", m, "U")


@needs_opensees
def test_a_crash_of_opensees_cannot_kill_the_caller():
    """A vecxz parallel to the member axis terminates the OpenSees interpreter (measured).  Through the
    runner the caller gets a failed result instead."""
    job = {"solver": "opensees", "secondOrder": False, "tolerance": 1e-10, "maxIterations": 10, "system": "UmfPack",
           "model": [["node", 1, 0.0, 0.0, 0.0], ["node", 2, 1000.0, 0.0, 0.0], ["fix", 1, 1, 1, 1, 1, 1, 1],
                     ["geomTransf", "Linear", 1, 1.0, 0.0, 0.0],
                     ["element", "elasticBeamColumn", 1, 1, 2, 5e3, 2e5, 7.7e4, 3e5, 8e6, 4e7, 1]],
           "pattern": [["load", 2, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]], "nodes": [1, 2], "elements": [1], "beams": [1]}
    out = run_job(job)
    assert out["ok"] is False and out["error"]
    assert "returncode" in out or "OpenSeesError" in out["error"] or "stderr" in out["error"]


@needs_opensees
def test_timeout_and_failed_analysis_are_results():
    m = _model()
    job = ops_tr.translate(m, "U").job
    out = run_job(job, timeout=0.001)
    assert out["ok"] is False and out.get("timeout") and "killed" in out["error"]
    # a mechanism (no supports at all in the foreign model) -> analyze fails, reported with the reason
    bad = dict(job, model=[c for c in job["model"] if c[0] != "fix"], system="FullGeneral")
    out = run_job(bad)
    assert out["ok"] is False and "singular" in out["error"]
    t = ops_tr.translate(m, "U")
    res = t.results(out)
    assert not res.ok and res.error and res.displacements is None
    # the default sparse solver (UmfPack) does NOT report the singular matrix (measured); the
    # translator's own self-check (finite, sane displacements; joint equilibrium) turns it into a failure
    out = run_job(dict(bad, system="UmfPack"))
    res = t.results(out)
    assert not res.ok and res.displacements is None
    assert "singular" in res.error or "self-check failed" in res.error


@needs_pynite
def test_worker_command_line(tmp_path):
    import subprocess

    m = _model()
    t = pyn_tr.translate(m, "U")
    src, dst = tmp_path / "job.json", tmp_path / "out.json"
    src.write_text(json.dumps(t.job))
    env = dict(__import__("os").environ)
    env["PYTHONPATH"] = str(importlib.import_module("benchfea").__path__[0] + "/..")
    proc = subprocess.run([sys.executable, "-m", "benchfea.crosscheck.worker", str(src), str(dst)], env=env,
                          capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0
    res = t.results(json.loads(dst.read_text()))
    assert res.ok and res.fidelity.level == EXACT and math.isfinite(float(res.displacements.sum()))
    assert equilibrium(m, res)["force"] < 1e-10
    assert subprocess.run([sys.executable, "-m", "benchfea.crosscheck.worker"], env=env, capture_output=True).returncode == 2
