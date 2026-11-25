"""Splits and the dataset builder (benchfea.bench.splits, benchfea.bench.dataset), the ``benchfea score``
and ``benchfea bench`` commands, and the dataset committed under data/.  docs/BENCHMARK.md sections 7-8."""
import json
from pathlib import Path

import pytest
from bench_helpers import CACHE, CELLS

from benchfea.bench import curveballs, dataset, splits
from benchfea.bench import scoring as S
from benchfea.bench import tasks as T
from benchfea.bench.families import FAMILIES
from benchfea.cli import main

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


# ------------------------------------------------------------------------------------------ splits
def test_split_file_is_well_formed():
    d = splits.load_splits()
    assert d["format"] == "benchfea-splits" and d["seedRule"] == {"modulus": 5, "heldOutResidues": [4]}
    held = [tuple(p) for p in d["heldOutCombinations"]]
    assert len(set(held)) == len(held) >= 8
    assert {f for f, _ in held} == set(FAMILIES)                  # every family has held-out combinations
    for f, c in held + [tuple(p) for p in d["servedCombinations"]]:
        assert f in FAMILIES and c in curveballs.available()
    served = {tuple(p) for p in d["servedCombinations"]}
    assert len(served) == len(d["servedCombinations"])
    assert served & set(held), "no held-out combination is served"
    assert served - set(held), "no train combination is served"


def test_is_held_out_rule():
    d = splits.load_splits()
    hf, hc = d["heldOutCombinations"][0]
    assert all(splits.is_held_out(hf, hc, s) for s in range(20))              # a held-out pair: every seed
    other = next(c for c in curveballs.available() if (hf, c) not in d["_heldOut"])
    assert [s for s in range(20) if splits.is_held_out(hf, other, s)] == [4, 9, 14, 19]   # else: seed % 5 == 4
    assert splits.split_of(hf, other, 3) == "train" and splits.split_of(hf, other, 4) == "held-out"
    assert splits.split_of(hf, hc, 3) == "held-out"
    assert splits.seeds("train", hf, other, 6) == [0, 1, 2, 3, 5, 6]
    assert splits.seeds("held-out", hf, other, 3) == [4, 9, 14]
    assert splits.seeds("train", hf, hc, 3) == [] and splits.seeds("held-out", hf, hc, 3) == [0, 1, 2]
    train = splits.combinations("train")
    assert train and not [p for p in train if tuple(p) in d["_heldOut"]]
    assert set(splits.combinations("held-out")) == set(splits.combinations()) > set(train)
    with pytest.raises(ValueError):
        splits.combinations("test")


# ----------------------------------------------------------------------------------------- builder
def test_plan_is_deterministic_and_covers_every_level():
    cells = dataset.plan(seeds=(0, 1, 4))
    assert cells == dataset.plan(seeds=(0, 1, 4))
    assert len(cells) == len(FAMILIES) * len(curveballs.available()) * 3 == len(set(cells))
    for family in FAMILIES:
        assert {c[3] for c in cells if c[0] == family} == set(FAMILIES[family].LEVELS)
    assert all(level in FAMILIES[f].LEVELS for f, _, _, level in cells)
    assert dataset.plan(seeds=(0,), families=["stair_frame"], curve_balls=["load_moved"]) == [
        ("stair_frame", "load_moved", 0, dataset.level_for("stair_frame", "load_moved", 0))]
    with pytest.raises(ValueError):
        dataset.plan(families=["castle"])
    with pytest.raises(KeyError):
        dataset.plan(curve_balls=["nope"])


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("dataset")
    cells = [CELLS[k] for k in ("change", "mechanism", "design_only", "held_out_seed")]
    manifest = dataset.build(root, cells, cache_dir=CACHE)
    return root, cells, manifest


def test_build_writes_only_verified_tasks_and_logs_the_rest(built):
    root, cells, m = built
    assert m["format"] == "benchfea-dataset" and "It does not mean reviewed" in m["meaning"]
    c = m["counts"]
    assert c["cells"] == 4 and c["cellsVerified"] == 3
    # change: predict, diagnose, repair; the two mechanism cells: diagnose, repair; design-only: nothing
    assert c["tasks"] == 7 == len(m["tasks"]) and c["byTaskType"] == {"diagnose": 3, "predict": 1, "repair": 3}
    assert c["bySplit"] == {"held-out": 2, "train": 5} and "design" not in c["byTaskType"]
    assert sum(c["byMode"].values()) == 7 and set(c["byMode"]) <= set(T.MODES)
    for t in m["tasks"]:
        task_file = root / "tasks" / t["file"]
        sub = dataset.SPLIT_DIRS[t["split"]]
        assert task_file.parent.name == sub and task_file.exists()
        public = json.loads(task_file.read_text())
        assert "model" not in public and (task_file.parent / public["modelFile"]).exists()
        assert not {"changeRecord", "reference", "key"} & set(public)
        task = dataset.load_task(task_file)
        answers = dataset.load_answers(task_file)
        assert (root / "answers" / sub / f"{t['id']}.answers.json").exists()
        assert answers["id"] == task["id"] == t["id"] and answers["verification"]["status"] == "verified"
        assert answers["changeRecord"]["curveBall"] and task["mode"] == dataset.mode_for(answers["key"], t["taskType"])
        assert splits.split_of(answers["key"]["family"], answers["key"]["curveBall"], answers["key"]["seed"]) == t["split"]
        assert S.score_answer(task, answers, S.reference_answer(task, answers))["score"] == 1.0
    assert [p.name for p in dataset.iter_tasks(root, "held-out")] == sorted(
        t["file"].split("/")[1] for t in m["tasks"] if t["split"] == "held-out")
    # nothing of the answers under tasks/
    for f in (root / "tasks").rglob("*.json"):
        assert "changeRecord" not in f.read_text() or f.name == "manifest.json"
    assert "changeRecord" not in (root / "tasks" / "manifest.json").read_text()
    # the quarantine log: the design-only cell (every attempt), and the task types not served on the others
    log = m["quarantineLog"]
    attempts = [e for e in log if e["what"] == "attempt"]
    assert len(attempts) == T.MAX_ATTEMPTS and {e["curveBall"] for e in attempts} == {CELLS["design_only"][1]}
    assert all(e["outcome"] == "cell not served" and e["reason"] for e in attempts)
    assert any(e["category"] == "not-covered" for e in attempts)
    not_served = {(e["curveBall"], e["taskType"]) for e in log if e["what"] == "task type"}
    assert not_served == {("release_added", "design"), ("release_mechanism", "predict"), ("release_mechanism", "design")}
    pair = m["cellsByFamilyAndCurveBall"]["continuous_beam|design_bracing_points_changed"]
    assert pair["verified"] == 0 and pair["cells"] == 1 and sum(pair["reasons"].values()) == 1
    assert m["cellsByFamily"]["continuous_beam"] == {"cells": 4, "verified": 3, "tasks": 7,
                                                     "notVerified": pair["reasons"]}
    assert m["cellsByCurveBall"]["release_mechanism"]["verified"] == 2


def test_build_is_deterministic_resumable_and_parallel(built, tmp_path):
    root, cells, m = built
    again = dataset.build(tmp_path, cells, cache_dir=CACHE, workers=2)        # from the cache, two processes

    def stable(manifest):
        return json.dumps({k: v for k, v in manifest.items() if k != "versions"}, sort_keys=True)

    assert stable(again) == stable(m)
    for f in sorted((root / "tasks").rglob("*.json")) + sorted((root / "answers").rglob("*.json")):
        assert (tmp_path / f.relative_to(root)).read_text() == f.read_text()
    # clean: files this run did not write are removed
    stray = tmp_path / "tasks" / "train" / "t00000000000000.task.json"
    stray.write_text("{}")
    dataset.build(tmp_path, cells[:1], cache_dir=CACHE, clean=True)
    assert not stray.exists() and len(list((tmp_path / "tasks" / "train").glob("*.task.json"))) == 3
    assert not list((tmp_path / "tasks" / "heldout").glob("*.json"))


def test_environment_can_serve_the_stored_dataset(built):
    """BenchFEAEnv(dataset_dir=...) serves the stored tasks -- the same task the environment builds on the
    fly for the same instance -- without generating or verifying anything at reset."""
    pytest.importorskip("gymnasium")
    from benchfea.bench.env import BenchFEAEnv, run_episode
    from benchfea.bench.tasks import TaskUnavailable
    root, cells, m = built
    stored, live = BenchFEAEnv(dataset_dir=root), BenchFEAEnv(cache_dir=CACHE)
    seen = set()
    for seed in range(8):
        obs, info = stored.reset(seed=seed)
        assert info["split"] == "train" and info["task_id"] in {t["id"] for t in m["tasks"] if t["split"] == "train"}
        seen.add(info["task_id"])
        obs2, info2 = stored.reset(seed=seed)
        assert obs == obs2 and info == info2                       # deterministic per seed
        key = stored.hidden_answers()["key"]
        o, i = live.reset(seed=0, options={"family": key["family"], "curve_ball": key["curveBall"],
                                           "instance_seed": key["seed"], "difficulty": key["difficulty"],
                                           "task_type": info["task_type"], "mode": info["mode"]})
        assert o == obs and i == info and live.hidden_answers() == stored.hidden_answers()
    assert len(seen) >= 3
    obs, info = stored.reset(seed=0, options={"task_type": "repair", "curve_ball": "release_mechanism"})
    assert info["task_type"] == "repair"
    held = BenchFEAEnv(dataset_dir=root, split="held-out")
    obs, info = held.reset(seed=0)
    assert info["split"] == "held-out" and info["task_id"] in {t["id"] for t in m["tasks"] if t["split"] == "held-out"}
    with pytest.raises(TaskUnavailable):
        stored.reset(seed=0, options={"family": "pipe_rack"})
    with pytest.raises(ValueError, match="held-out"):
        stored.reset(seed=0, options={"family": "continuous_beam", "curve_ball": "release_mechanism", "instance_seed": 4})

    def oracle(observation, info):
        from benchfea.bench.scoring import reference_answer
        answer = reference_answer({"taskType": info["task_type"]}, stored.hidden_answers())
        return json.dumps({"tool": "submit", "args": {"answer": answer}})

    assert run_episode(stored, oracle, seed=3)["reward"] == 1.0


# --------------------------------------------------------------------------------------------- CLI
def test_cli_score(built, tmp_path, capsys):
    root, cells, m = built
    for t in m["tasks"][:3]:
        task_file = root / "tasks" / t["file"]
        task, answers = dataset.load_task(task_file), dataset.load_answers(task_file)
        good, bad = tmp_path / "good.json", tmp_path / "bad.json"
        good.write_text(json.dumps(S.reference_answer(task, answers)))
        bad.write_text("this is not json")
        assert main(["score", str(good), "--task", str(task_file), "--json"]) == 0
        out = json.loads(capsys.readouterr().out)
        assert out["score"] == 1.0 and out["taskType"] == t["taskType"] and out["format"] == "benchfea-score"
        assert main(["score", str(bad), "--task", str(task_file), "-o", str(tmp_path / "s.json")]) == 0
        text = capsys.readouterr().out
        assert "score 0.0000 (answer not valid)" in text and "not valid JSON" in text
        assert json.loads((tmp_path / "s.json").read_text())["score"] == 0.0
        # explicit answers file; text output
        answers_file = root / "answers" / task_file.parent.name / f"{t['id']}.answers.json"
        assert main(["score", str(good), "--task", str(task_file), "--answers", str(answers_file)]) == 0
        assert "score 1.0000" in capsys.readouterr().out
    other = root / "answers" / "train" / f"{m['tasks'][-1]['id']}.answers.json"
    other = other if other.exists() else root / "answers" / "heldout" / f"{m['tasks'][-1]['id']}.answers.json"
    assert main(["score", str(good), "--task", str(root / "tasks" / m["tasks"][0]["file"]), "--answers", str(other)]) == 2
    assert main(["score", str(good), "--task", str(tmp_path / "missing.task.json")]) == 2
    assert main(["score", str(tmp_path / "missing.json"), "--task", str(root / "tasks" / m["tasks"][0]["file"])]) == 2
    capsys.readouterr()


def test_cli_bench_make_dataset_and_run_agent(tmp_path, capsys):
    assert main(["bench", "make-dataset", "--out", str(tmp_path), "--seeds", "0", "--families", "continuous_beam",
                 "--curve-balls", "release_mechanism", "--cache", CACHE, "--quiet"]) == 0
    assert "2 verified task(s) from 1 of 1 cell(s)" in capsys.readouterr().out
    assert json.loads((tmp_path / "tasks" / "manifest.json").read_text())["counts"]["tasks"] == 2
    assert main(["bench", "make-dataset", "--families", "castle"]) == 2
    assert main(["bench"]) == 2
    pytest.importorskip("gymnasium")
    agents = ROOT / "examples" / "agents"
    common = ["--family", "continuous_beam", "--curve-ball", "release_mechanism", "--difficulty", "1", "--cache", CACHE,
              "--json"]
    capsys.readouterr()
    assert main(["bench", "run-agent", str(agents / "oracle_agent.py"), "--episodes", "2", "--seed", "0"] + common) == 0
    out = json.loads(capsys.readouterr().out)
    assert out["meanReward"] == 1.0 and len(out["episodes"]) == 2 and all(e["terminated"] for e in out["episodes"])
    assert main(["bench", "run-agent", str(agents / "random_agent.py")] + common) == 0
    assert json.loads(capsys.readouterr().out)["meanReward"] < 0.5
    assert main(["bench", "run-agent", str(agents / "missing.py")]) == 2
    assert main(["bench", "run-agent", str(agents / "oracle_agent.py"), "--task-type", "predict"] + common) == 8
    capsys.readouterr()


# ------------------------------------------------------------------------------- committed dataset
@pytest.fixture(scope="module")
def committed():
    manifest = DATA / "tasks" / "manifest.json"
    assert manifest.exists(), "the dataset under data/ is missing: benchfea bench make-dataset"
    return json.loads(manifest.read_text())


def test_committed_dataset_matches_its_manifest(committed):
    m = committed
    assert m["counts"]["tasks"] == len(m["tasks"]) >= 200
    listed = {t["id"]: t for t in m["tasks"]}
    assert len(listed) == len(m["tasks"])
    on_disk = {}
    for split, sub in dataset.SPLIT_DIRS.items():
        for f in dataset.iter_tasks(DATA, split):
            on_disk[f.name.split(".")[0]] = (split, f)
        models = {f.name.split(".")[0] for f in (DATA / "tasks" / sub).glob("*.model.json")}
        answers = {f.name.split(".")[0] for f in (DATA / "answers" / sub).glob("*.answers.json")}
        tasks = {k for k, (s, _) in on_disk.items() if s == split}
        assert models == tasks == answers
    assert set(on_disk) == set(listed)
    held = splits.load_splits()["_heldOut"]
    for tid, (split, f) in on_disk.items():
        public = json.loads(f.read_text())
        assert public["split"] == split == listed[tid]["split"] and public["taskType"] == listed[tid]["taskType"]
        assert not {"changeRecord", "reference", "key"} & set(public)
        assert public["verification"]["status"] == "verified"
    for split, sub in dataset.SPLIT_DIRS.items():
        for f in sorted((DATA / "answers" / sub).glob("*.answers.json")):
            a = json.loads(f.read_text())
            key = a["key"]
            assert splits.split_of(key["family"], key["curveBall"], key["seed"]) == split
            assert (split == "held-out") or (key["family"], key["curveBall"]) not in held
            assert a["verification"]["status"] == "verified" and a["verification"]["tiersPassed"]
            assert not a["verification"]["tiersNotPassed"]
    # both splits, every family, the three task types that can be served, both modes, several levels
    c = m["counts"]
    assert set(c["bySplit"]) == {"train", "held-out"} and set(c["byFamily"]) == set(FAMILIES)
    assert {"predict", "diagnose", "repair"} <= set(c["byTaskType"]) and set(c["byMode"]) == set(T.MODES)
    assert len(c["byDifficulty"]) >= 3
    # the served combinations of the split file are those with a verified cell in the manifest
    served = {tuple(k.split("|")) for k, v in m["cellsByFamilyAndCurveBall"].items() if v["verified"]}
    assert served == {tuple(p) for p in splits.load_splits()["servedCombinations"]}
    # some combinations appear only in the held-out split
    train_pairs = {(json.loads((DATA / "answers" / "train" / f"{t['id']}.answers.json").read_text())["key"]["family"],
                    json.loads((DATA / "answers" / "train" / f"{t['id']}.answers.json").read_text())["key"]["curveBall"])
                   for t in m["tasks"] if t["split"] == "train"}
    assert (served & held) and not (train_pairs & held)


def test_committed_tasks_are_solved_by_their_reference_answers(committed):
    """Every predict and diagnose task of the committed dataset scores 1.0 with the answer built from its
    hidden reference; a sample of the repair tasks is re-solved through the verification pipeline."""
    repairs = 0
    for t in committed["tasks"]:
        f = DATA / "tasks" / t["file"]
        if t["taskType"] == "repair":
            if repairs >= 2 or t["family"] != "continuous_beam":
                continue
            repairs += 1
        task, answers = dataset.load_task(f), dataset.load_answers(f)
        assert S.score_answer(task, answers, S.reference_answer(task, answers))["score"] == 1.0, t["id"]
    assert repairs == 2


def test_committed_task_is_reproduced_from_its_key(committed):
    """One committed task rebuilt from nothing: same model, same reference, same verification outcome."""
    t = next(t for t in committed["tasks"] if t["family"] == "continuous_beam" and t["taskType"] == "predict"
             and t["difficulty"] == 1)
    f = DATA / "tasks" / t["file"]
    task, answers = dataset.load_task(f), dataset.load_answers(f)
    key = answers["key"]
    case = T.build_case(key["family"], key["curveBall"], key["seed"], key["difficulty"])
    rebuilt, rebuilt_answers = T.make_task(case, task["taskType"], task["mode"])
    task.pop("modelFile")                                         # only the stored file has it
    assert rebuilt == task
    assert rebuilt_answers["reference"] == answers["reference"]
    assert rebuilt_answers["verification"]["tiersPassed"] == answers["verification"]["tiersPassed"]
    assert rebuilt_answers["changeRecord"]["changes"] == answers["changeRecord"]["changes"]
