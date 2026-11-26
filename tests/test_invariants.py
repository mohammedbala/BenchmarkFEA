"""Invariants reported on every solve, the linearity checks (superposition, Maxwell-Betti), tests that
the checks have teeth (injected faults are caught), determinism, and the separation of design-only data
from the solver."""
import json
import pathlib
import re

import numpy as np
import pytest

import benchfea.solve
from benchfea.model import build_model
from benchfea.solve import (LinearAnalysis, ModelRefused, analyze, reciprocity_check, superposition_check)
from benchfea.solve.element import MemberElement
from benchfea.solve.invariants import rigid_body_modes

from solve_reference import make

NAMES = ["stiffness symmetry", "rigid-body modes", "positive definiteness", "linear solve residual",
         "global force equilibrium", "global moment equilibrium", "joint equilibrium", "member equilibrium",
         "member compatibility", "work = strain energy"]


def space_frame(design=None, extra_cases=None):
    """Two-bay, one-storey space frame using every feature of the linear core: rolled and skew members,
    full and partial releases, rigid offsets (also off-axis), a rigid link, a diaphragm, springs, a truss
    brace, every member load type, self weight, a prescribed support displacement."""
    H = 3500.0
    nodes = {"a1": (0, 0, 0), "a2": (6000, 0, 0), "a3": (6000, 5000, 0), "a4": (0, 5000, 0), "a5": (11000, 800, 0),
             "t1": (0, 0, H), "t2": (6000, 0, H), "t3": (6000, 5000, H), "t4": (0, 5000, H), "t5": (11000, 800, H + 900),
             "cm": (3000, 2500, H), "ecc": (6400, 150, H), "mid": (3000, 0, H)}
    members = {
        "c1": {"i": "a1", "j": "t1", "orientation": {"roll": 30}},
        "c2": {"i": "a2", "j": "t2"},
        "c3": {"i": "a3", "j": "t3", "orientation": {"roll": 90}},
        "c4": {"i": "a4", "j": "t4", "releases": {"i": {"ry": {"spring": 4.0e9}}}},
        "c5": {"i": "a5", "j": "t5", "orientation": {"referenceVector": [0, 1, 0]}},
        "g1a": {"i": "t1", "j": "mid", "offsets": {"i": [150, 0, -120]}},
        "g1b": {"i": "mid", "j": "t2", "releases": {"j": {"rz": True}}},
        "g2": {"i": "t2", "j": "t3", "releases": {"i": {"rz": {"spring": 6.0e9}}, "j": {"rz": {"spring": 1.0e7}}}},
        "g3": {"i": "t4", "j": "t3"},
        "g4": {"i": "t1", "j": "t4", "orientation": {"roll": 10}},
        "rafter": {"i": "ecc", "j": "t5", "offsets": {"j": [-200, 0, -100]}, "releases": {"j": {"ry": True, "ux": {"spring": 5.0e4}}}},
        "brace": {"i": "a2", "j": "t5", "type": "truss"},
    }
    supports = {"a1": "fixed", "a2": "pinned", "a3": "fixed", "a4": {"type": "pinned", "dofs": {"rz": "fixed"}},
                "a5": {"type": "fixed", "dofs": {"ux": {"spring": 9000.0}, "ry": {"spring": 5.0e9}}}}
    constraints = {"floor": {"type": "rigidDiaphragm", "master": "cm", "slaves": ["t1", "t2", "t3", "t4", "mid"]},
                   "bracket": {"type": "rigidLink", "master": "ecc", "slaves": []}}
    # the bracket: node "ecc" is rigidly attached to the top of column c2 ... t2 is a diaphragm slave already
    # (a node may be the slave of one constraint only), so the link runs from ecc (master) to a helper node
    nodes["helper"] = (6400, 150, H + 500)
    members["post"] = {"i": "helper", "j": "t5"}
    constraints["bracket"]["slaves"] = ["helper"]
    members["stub"] = {"i": "t2", "j": "ecc"}
    cases = {
        "D": [{"type": "selfWeight"},
              {"type": "memberDistributed", "member": "g1a", "direction": "Z", "w1": -8.0},
              {"type": "memberDistributed", "member": "g1b", "direction": "Z", "w1": -8.0, "w2": -3.0, "start": 0.1, "end": 0.8},
              {"type": "memberDistributed", "member": "rafter", "direction": "Z", "w1": -4.0, "projected": True}],
        "L": [{"type": "memberPoint", "member": "g2", "direction": "Z", "value": -20000.0, "position": 0.4},
              {"type": "memberPoint", "member": "g3", "direction": "y", "kind": "moment", "value": 6.0e6, "position": 0.7},
              {"type": "memberPoint", "member": "c1", "direction": "x", "kind": "moment", "value": 2.0e6, "position": 0.5},
              {"type": "memberDistributed", "member": "brace", "direction": "Y", "w1": 0.5},
              {"type": "nodal", "node": "helper", "Fx": 3000.0, "Mz": 4.0e6}],
        "W": [{"type": "nodal", "node": "cm", "Fx": 40000.0, "Fy": -25000.0, "Mz": 5.0e7},
              {"type": "memberDistributed", "member": "c1", "direction": "x", "w1": 1.0},
              {"type": "memberDistributed", "member": "c5", "direction": "X", "w1": 2.0, "w2": 0.0}],
        "S": [{"type": "supportDisplacement", "node": "a3", "uz": -8.0, "rx": 0.002}],
    }
    cases.update(extra_cases or {})
    combos = {"U1": {"factors": {"D": 1.4}}, "U2": {"factors": {"D": 1.2, "L": 1.6}},
              "U3": {"factors": {"D": 1.2, "L": 0.5, "W": 1.0, "S": 1.0}}, "U4": {"factors": {"D": 0.9, "W": -1.0}},
              "SLS": {"factors": {"D": 1.0, "L": 1.0, "S": 1.0}, "purpose": "serviceability"}}
    for k in (extra_cases or {}):
        combos[k] = {"factors": {k: 1.0}}
    return make(nodes, members, supports, cases, combos=combos, constraints=constraints, design=design,
                materials={"M": {"E": 200000.0, "G": 77000.0, "unitWeight": 7.7e-5, "Fy": 345.0, "Fu": 450.0}})


@pytest.fixture(scope="module")
def frame():
    model = build_model(space_frame())
    return model, analyze(model)


def test_every_combination_carries_the_full_invariants_report_and_all_pass(frame):
    model, res = frame
    tol = model.analysis.equilibrium_tolerance
    for c in res.combinations.values():
        assert [ch.name for ch in c.invariants.checks] == NAMES
        assert c.invariants.passed
        by = {ch.name: ch for ch in c.invariants.checks}
        assert by["stiffness symmetry"].tolerance == 1e-10 and by["rigid-body modes"].tolerance == 1e-9
        assert by["positive definiteness"].tolerance == model.analysis.pivot_tolerance
        assert by["positive definiteness"].value > by["positive definiteness"].tolerance
        for n in NAMES[4:]:
            assert by[n].tolerance == tol and 0.0 <= by[n].value < 1e-9      # far inside the tolerance
        assert by["linear solve residual"].value < 1e-12
        assert all(ch.detail for ch in c.invariants.checks)


def test_global_equilibrium_recomputed_here_from_the_file_loads(frame):
    """Applied loads summed in the test straight from the model's load objects (points and line-load
    resultants by numerical quadrature), then added to the reported reactions."""
    from scipy.integrate import quad
    from benchfea.model import MemberDistributedLoad, MemberPointLoad, NodalLoad, SelfWeight
    model, res = frame

    def case_resultant(case):
        F, M = np.zeros(3), np.zeros(3)
        for ld in case.loads:
            if isinstance(ld, NodalLoad):
                F += ld.values[:3]
                M += ld.values[3:] + np.cross(model.nodes[ld.node].xyz, ld.values[:3])
            elif isinstance(ld, MemberPointLoad):
                m = model.members[ld.member]
                v = m.R.T @ ld.vector
                if ld.kind == "force":
                    F += v
                    M += np.cross(m.start + ld.x * m.R[0], v)
                else:
                    M += v
            elif isinstance(ld, (MemberDistributedLoad, SelfWeight)):
                items = []
                if isinstance(ld, SelfWeight):
                    for mid in ld.members:
                        m = model.members[mid]
                        w = ld.factor * model.materials[m.material].unit_weight * model.sections[m.section].A * ld.direction
                        items.append((m, 0.0, m.length, w, w))
                else:
                    m = model.members[ld.member]
                    items.append((m, ld.x1, ld.x2, m.R.T @ ld.w1, m.R.T @ ld.w2))
                for m, x1, x2, w1, w2 in items:
                    w = lambda x: w1 + (w2 - w1) * (x - x1) / (x2 - x1)
                    for c in range(3):
                        F[c] += quad(lambda x: w(x)[c], x1, x2, epsabs=0, epsrel=1e-13)[0]
                        M[c] += quad(lambda x: np.cross(m.start + x * m.R[0], w(x))[c], x1, x2, epsabs=0, epsrel=1e-13)[0]
        return F, M

    per_case = {cid: case_resultant(case) for cid, case in model.load_cases.items()}
    xyz = np.array([n.xyz for n in model.nodes.values()])
    for cid, combo in model.combinations.items():
        F = sum(f * per_case[k][0] for k, f in combo.factors.items())
        M = sum(f * per_case[k][1] for k, f in combo.factors.items())
        R = res.combinations[cid].reactions
        Fr = R[:, :3].sum(axis=0)
        Mr = (R[:, 3:] + np.cross(xyz, R[:, :3])).sum(axis=0)
        assert np.abs(F + Fr).max() <= 1e-9 * np.abs(R[:, :3]).sum()
        assert np.abs(M + Mr).max() <= 1e-9 * (np.abs(R[:, 3:]) + np.abs(np.cross(xyz, R[:, :3]))).sum()


def test_joint_equilibrium_recomputed_here_from_member_end_forces(frame):
    """At every node: nodal loads + reactions + constraint forces = sum of the member end forces carried
    to the node (rotation to global and lever arm of the rigid offset), with the transformation written
    out here."""
    from benchfea.model import NodalLoad
    model, res = frame
    for cid, combo in model.combinations.items():
        c = res.combinations[cid]
        g = c.reactions + c.constraint_forces
        for case_id, f in combo.factors.items():
            for ld in model.load_cases[case_id].loads:
                if isinstance(ld, NodalLoad):
                    g[model.nodes[ld.node].index] += f * ld.values
        scale_f = scale_m = 0.0
        for mid, m in model.members.items():
            ef = c.members[mid].end_forces
            for blk, node, off in ((0, m.i, m.offset_i), (6, m.j, m.offset_j)):
                F = m.R.T @ ef[blk:blk + 3]
                M = m.R.T @ ef[blk + 3:blk + 6] + np.cross(off, F)        # force at the flexible end, arm = offset
                g[model.nodes[node].index, :3] -= F
                g[model.nodes[node].index, 3:] -= M
                scale_f, scale_m = max(scale_f, np.abs(F).max()), max(scale_m, np.abs(M).max())
        assert np.abs(g[:, :3]).max() <= 1e-9 * scale_f and np.abs(g[:, 3:]).max() <= 1e-9 * scale_m


def test_stiffness_matrix_properties_checked_directly(frame):
    model, _ = frame
    la = LinearAnalysis(model)
    K = la.K_members
    assert np.abs(K - K.T).max() <= 1e-12 * np.abs(K).max()
    B = rigid_body_modes(la.structure.xyz)
    # zero strain energy and zero forces under the six rigid-body motions of the unsupported structure
    energy = np.einsum("ij,ij->j", B, K @ B)
    assert np.abs(energy).max() <= 1e-10 * np.abs(np.einsum("ij,ij->j", np.abs(B), np.abs(K) @ np.abs(B))).max()
    lam = np.linalg.eigvalsh(0.5 * (K + K.T))
    assert lam.min() >= -1e-9 * lam.max()                                  # positive semi-definite
    # after supports and constraints: positive definite (Cholesky succeeds on the reduced matrix)
    st = la.structure
    Kr = st.reduce(la.K)
    Kff = Kr[np.ix_(st.free_r, st.free_r)]
    np.linalg.cholesky(Kff)
    f = la.factor
    assert f.min_pivot > 0 and f.lambda_min > 0 and f.condition == pytest.approx(f.lambda_max / f.lambda_min)
    # constraints are exact: T reproduces the rigid relations, and slaves are not unknowns
    assert Kff.shape[0] == st.n - st.slave.sum() - st.fixed.sum() - st.inactive.sum()


def test_superposition_holds(frame):
    model, _ = frame
    checks = superposition_check(model)
    assert [c.name for c in checks] == [f"superposition: {k}" for k in model.combinations]
    assert all(c.passed and c.value < 1e-10 for c in checks)


def test_reciprocity_maxwell_betti(frame):
    model, _ = frame
    auto = reciprocity_check(model)
    assert auto.passed and auto.value < 1e-10 and "12 DOF pairs" in auto.detail
    pairs = [(("t5", "uz"), ("mid", "uz")), (("t1", "rx"), ("helper", "uy")), (("t3", "ux"), ("t5", "ry")),
             (("cm", "rz"), ("t5", "ux")), (("t2", "uy"), ("ecc", "rz"))]       # includes slave and master DOFs
    chosen = reciprocity_check(model, pairs)
    assert chosen.passed and chosen.value < 1e-10


def test_reciprocity_through_the_public_analysis_of_two_load_cases():
    """Unit force at A (direction p) and unit moment at B (direction q), each as its own load case:
    rotation at B about q due to the force = displacement at A along p due to the moment."""
    extra = {"PA": [{"type": "nodal", "node": "t5", "Fz": 1.0}], "MB": [{"type": "nodal", "node": "t4", "My": 1.0}]}
    model = build_model(space_frame(extra_cases=extra))
    res = analyze(model)
    a, b = model.nodes["t5"].index, model.nodes["t4"].index
    d_ba = res.combinations["PA"].displacements[b, 4]
    d_ab = res.combinations["MB"].displacements[a, 2]
    assert d_ba != 0.0 and d_ab == pytest.approx(d_ba, rel=1e-9)


def test_results_are_deterministic_bit_for_bit(frame):
    model, res = frame
    again = analyze(build_model(space_frame()))
    assert json.dumps(res.to_document(model)) == json.dumps(again.to_document(model))


# ------------------------------------------------------------------ the checks have teeth
def simple():
    return make({"a": (0, 0, 0), "b": (4000, 0, 0), "c": (4000, 3000, 0)}, {"m": ("a", "b"), "n": ("b", "c")}, {"a": "fixed"},
                {"c": [{"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -5.0},
                       {"type": "nodal", "node": "c", "Fz": -2000.0}]})


def failed_names(exc):
    return sorted({re.search(r"invariant '([^']+)'", p.message).group(1) for p in exc.problems})


def test_wrong_fixed_end_forces_are_caught(monkeypatch):
    original = MemberElement._fixed_end
    monkeypatch.setattr(MemberElement, "_fixed_end", lambda self, loads: tuple(1.001 * x for x in original(self, loads)))
    with pytest.raises(ModelRefused) as e:
        analyze(build_model(simple()))
    assert e.value.codes and set(e.value.codes) == {"solve.invariant-failed"}
    assert "global force equilibrium" in failed_names(e.value) and "member equilibrium" in failed_names(e.value)
    flagged = analyze(build_model(simple()), strict=False)
    assert not flagged.combinations["c"].invariants.passed and flagged.combinations["c"].warnings
    assert any("FLAGGED" in w for w in flagged.warnings)
    assert e.value.results.combinations["c"].invariants.passed is False


def test_wrong_element_stiffness_in_recovery_is_caught(monkeypatch):
    """Assemble with one stiffness, recover with another: joint equilibrium and compatibility fail."""
    original = MemberElement.recover

    def recover(self, u, loads, stations):
        return original(self, 1.0005 * u, loads, stations)
    monkeypatch.setattr(MemberElement, "recover", recover)
    with pytest.raises(ModelRefused) as e:
        analyze(build_model(simple()))
    assert "joint equilibrium" in failed_names(e.value)


def test_unsymmetric_stiffness_is_caught(monkeypatch):
    import benchfea.solve.analysis as analysis
    original = analysis.build_element

    def build(model, member, *a, **k):
        el = original(model, member, *a, **k)
        el.k_global = el.k_global.copy()
        el.k_global[0, 7] *= 1.0 + 1e-6
        el.k_global[1, 5] += 1e-6 * abs(el.k_global[1, 1]) * 1000.0
        return el
    monkeypatch.setattr(analysis, "build_element", build)
    res = analyze(build_model(simple()), strict=False)
    names = [ch.name for ch in res.combinations["c"].invariants.checks if not ch.passed]
    assert "stiffness symmetry" in names


def test_element_that_resists_rigid_body_motion_is_caught(monkeypatch):
    import benchfea.solve.analysis as analysis
    original = analysis.build_element

    def build(model, member, *a, **k):
        el = original(model, member, *a, **k)
        el.k_global = el.k_global.copy()
        el.k_global[2, 2] *= 1.0 + 1e-5                    # a spurious spring to ground
        return el
    monkeypatch.setattr(analysis, "build_element", build)
    res = analyze(build_model(simple()), strict=False)
    names = [ch.name for ch in res.combinations["c"].invariants.checks if not ch.passed]
    assert "rigid-body modes" in names


# ------------------------------------------------------- design-only data and the solver
def test_solver_sources_never_touch_design_only_data():
    """No module of benchfea.solve reads Model.design_members or imports the design layer."""
    root = pathlib.Path(benchfea.solve.__file__).parent
    files = sorted(root.glob("*.py"))
    assert len(files) >= 8
    for f in files:
        text = f.read_text(encoding="utf-8")
        code = "\n".join(line.split("#", 1)[0] for line in text.splitlines())
        code = re.sub(r'"""(?:.|\n)*?"""', "", code)                       # docstrings may mention it
        assert "design_members" not in code, f.name
        assert not re.search(r"\bdesignparams\b|\bDesignParams\b|\bdeflection_limits\b|benchfea\.design\b|from \.\.design\b|from \.\. import design", code), f.name
        assert not re.search(r"\.(Fy|Fu|Kx|Ky|Kz|Lb|Cb)\b", code), f.name


def test_design_members_are_never_accessed_at_run_time():
    class Poison(dict):
        def __getattribute__(self, name):
            raise AssertionError(f"the solver touched Model.design_members ({name})")
    model = build_model(space_frame())
    model.design_members = Poison()
    res = analyze(model)
    superposition_check(model)
    reciprocity_check(model)
    assert res.combinations


def test_changing_design_variables_gives_bitwise_identical_analysis_results(frame):
    model, res = frame
    design = {
        "defaults": {"method": "ASD", "E": 150000.0, "G": 50000.0, "Fy": 250.0, "Fu": 400.0, "Kx": 2.5, "Ky": 0.7, "Cb": 1.3,
                     "Lb": 1234.0, "U": 0.6, "role": "brace", "sideswayInhibited": True, "Cmx": 0.6,
                     "maxSlendernessCompression": 50, "deflectionLimits": [{"spanRatio": 100}]},
        "groups": {"cols": {"members": ["c1", "c2", "c3"], "Ly": 500.0, "Kz": 3.0, "bracing": {"torsion": [0.3]}}},
        "members": {"girder": {"members": ["g3"], "Kx": 0.5, "Fy": 500.0}},
    }
    model2 = build_model(space_frame(design=design))
    assert model2.design_members.keys() != model.design_members.keys()
    assert model2.design_members["c1"].E == 150000.0 and model2.materials["M"].E == 200000.0
    res2 = analyze(model2)
    for cid, c in res.combinations.items():
        c2 = res2.combinations[cid]
        assert np.array_equal(c.displacements, c2.displacements) and np.array_equal(c.reactions, c2.reactions)
        assert np.array_equal(c.constraint_forces, c2.constraint_forces)
        for mid, m in c.members.items():
            m2 = c2.members[mid]
            assert np.array_equal(m.end_forces, m2.end_forces) and np.array_equal(m.x, m2.x)
            assert np.array_equal(m.forces, m2.forces) and np.array_equal(m.displacements, m2.displacements)
        assert [(ch.name, ch.value) for ch in c.invariants.checks] == [(ch.name, ch.value) for ch in c2.invariants.checks]
    d1, d2 = res.to_document(model), res2.to_document(model2)
    for key in ("combinations", "envelopes", "nodalEnvelopes", "buckling"):
        assert json.dumps(d1[key]) == json.dumps(d2[key])
    assert d1["defaultsUsed"] != d2["defaultsUsed"]                        # only the report differs
