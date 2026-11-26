"""The Gymnasium environment (benchfea.bench.env): interface, tools, determinism, reward and budget,
no leak of the reference before the episode ends, held-out policy.  docs/BENCHMARK.md sections 5-7."""
import importlib.util
import json
from pathlib import Path

import pytest

gymnasium = pytest.importorskip("gymnasium")

from bench_helpers import CACHE, CELLS, case  # noqa: E402

from benchfea.bench import env as E  # noqa: E402
from benchfea.bench import scoring as S  # noqa: E402
from benchfea.bench import splits  # noqa: E402
from benchfea.bench.tasks import TaskUnavailable, make_task  # noqa: E402

AGENTS = Path(__file__).resolve().parents[1] / "examples" / "agents"


def options(name, **more):
    family, cb, seed, level = CELLS[name]
    return {"family": family, "curve_ball": cb, "instance_seed": seed, "difficulty": level, **more}


def make(**kw):
    kw.setdefault("cache_dir", CACHE)
    return E.BenchFEAEnv(**kw)


def call(tool, **args):
    return json.dumps({"tool": tool, "args": args})


def result(obs):
    return json.loads(obs["tool_result"])


def load_agent(name):
    spec = importlib.util.spec_from_file_location(name, AGENTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------------------- interface
def test_registered_and_spaces():
    assert E.ENV_ID == "BenchFEA-v0" and E.ENV_ID in gymnasium.registry
    env = gymnasium.make(E.ENV_ID, cache_dir=CACHE, max_steps=5)
    base = env.unwrapped
    assert isinstance(base, E.BenchFEAEnv) and isinstance(base, gymnasium.Env)
    assert isinstance(base.action_space, gymnasium.spaces.Text)
    assert isinstance(base.observation_space, gymnasium.spaces.Dict)
    assert set(base.observation_space.spaces) == {"statement", "input_file", "tool_result", "steps_remaining"}
    obs, info = env.reset(seed=1, options=options("change", task_type="predict", mode="tool-allowed"))
    assert base.observation_space.contains(obs)
    assert base.action_space.contains(call("inspect", path="/nodes"))
    model = json.loads(obs["input_file"])
    assert model == case("change")["perturbed"]["document"]
    assert obs["statement"].startswith("TASK (predict)") and obs["steps_remaining"] == 5
    obs, reward, terminated, truncated, info = env.step(call("inspect", path="/units"))
    assert base.observation_space.contains(obs) and isinstance(reward, float)
    assert terminated is False and truncated is False and isinstance(info, dict)
    env.close()
    env.close()


def test_check_env():
    """gymnasium.utils.env_checker.check_env passes on the environment (the random draw is restricted to
    four small cases so that the test is fast), the render check included.  Its only warning: no render
    fps is declared (the render mode is "ansi": text)."""
    import warnings

    from gymnasium.utils.env_checker import check_env
    env = gymnasium.make(E.ENV_ID, cache_dir=CACHE, max_steps=6, families=("continuous_beam",),
                         curve_balls=("support_fixed_to_pinned", "support_removed"), difficulties=(2,), seed_range=2)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        check_env(env.unwrapped)
    messages = [str(w.message) for w in caught]
    assert all("render fps" in m for m in messages), messages
    env.close()


def test_reset_without_options_draws_a_task_of_the_split():
    env = make(families=("continuous_beam",), curve_balls=("support_fixed_to_pinned", "support_removed"), difficulties=(2,),
               seed_range=2)
    seen = set()
    for seed in range(6):
        obs, info = env.reset(seed=seed)
        assert info["split"] == "train" and info["family"] == "continuous_beam" and info["difficulty"] == 2
        assert info["task_type"] in ("predict", "diagnose", "repair") and info["mode"] in ("tool-allowed", "no-solver")
        seen.add((info["task_type"], info["mode"]))
    assert len(seen) >= 3


def test_reset_options_are_checked():
    env = make()
    for bad in ({"famly": "portal_frame"}, {"family": "castle"}, {"curve_ball": "nope"}, {"task_type": "guess"},
                {"mode": "x"}, {"family": "continuous_beam", "difficulty": 5}, {"instance_seed": -1},
                {"instance_seed": "3"}):
        with pytest.raises(ValueError):
            env.reset(seed=0, options=bad)
    with pytest.raises(ValueError):
        E.BenchFEAEnv(split="test")
    with pytest.raises(RuntimeError):
        make().step(call("inspect"))


# ----------------------------------------------------------------------------------- determinism
SCRIPT = [call("inspect", path="/supports"), call("run_analysis"), call("read_results", path="/combinations"),
          "not json", call("edit", ops=[{"op": "remove", "path": "/members"}]), call("inspect"),
          call("run_analysis", design=True), call("read_results", document="design", path="/status"),
          call("submit", answer={"values": {"q1": 1.0}})]


def transcript(env, seed, opts):
    obs, info = env.reset(seed=seed, options=opts)
    out = [(json.dumps(obs, sort_keys=True), json.dumps(info, sort_keys=True))]
    for action in SCRIPT:
        obs, reward, terminated, truncated, info = env.step(action)
        out.append((json.dumps(obs, sort_keys=True), reward, terminated, truncated,
                    json.dumps(info, sort_keys=True, default=str)))
        if terminated or truncated:
            break
    return out


def test_two_environments_with_the_same_seed_give_identical_transcripts():
    restrict = dict(families=("continuous_beam",), curve_balls=("support_fixed_to_pinned", "support_removed"),
                    difficulties=(2,), seed_range=2)
    a = transcript(make(**restrict), 11, None)
    b = transcript(make(**restrict), 11, None)                    # a second environment, nothing shared in memory
    assert a == b and len(a) >= 2
    c = transcript(make(cache_dir=None), 5, options("mechanism"))  # without the disk cache: built from nothing
    d = transcript(make(), 5, options("mechanism"))
    assert c == d
    seeds = {transcript(make(**restrict), s, None)[0] for s in range(5)}
    assert len(seeds) > 1                                         # other seeds give other tasks


# ------------------------------------------------------------------------------------------- tools
def test_tools_and_reward_only_on_submit():
    env = make(max_steps=40)
    obs, info = env.reset(seed=0, options=options("change", task_type="predict", mode="tool-allowed"))
    given = json.loads(obs["input_file"])
    rewards = []

    def step(action):
        obs, reward, terminated, truncated, info = env.step(action)
        rewards.append(reward)
        assert not truncated
        return result(obs), terminated, info

    r, _, _ = step(call("inspect", path="/members"))
    assert r["ok"] and r["value"] == given["members"]
    r, _, _ = step(call("inspect"))
    assert r["ok"] and r["value"] == given
    r, _, _ = step(call("inspect", path="/members/none"))
    assert not r["ok"] and "none" in r["error"]
    r, _, _ = step(call("read_results", path="/combinations"))
    assert not r["ok"] and "run_analysis first" in r["error"]
    r, _, _ = step(call("run_analysis"))
    assert r["ok"] and set(r["combinations"]) == set(given["combinations"])
    r, _, _ = step(call("read_results", path="/unitLabels"))
    assert r["ok"] and set(r["value"]) == {"length", "rotation", "force", "moment"}
    values = {}
    for q in info["quantities"]:
        if q["kind"] == "memberMaxMoment":
            r, _, _ = step(call("read_results", path=f"/combinations/{q['combination']}/members/{q['member']}/stations/{q['component']}"))
            values[q["id"]] = max(abs(x) for x in r["value"])
        else:
            table = {"displacement": "displacements", "reaction": "reactions"}.get(q["kind"])
            path = (f"/combinations/{q['combination']}/{table}/{q['node']}/{q['component']}" if table else
                    f"/combinations/{q['combination']}/members/{q['member']}/endForces/{q['end']}/{q['component']}")
            r, _, _ = step(call("read_results", path=path))
            values[q["id"]] = r["value"]
    r, _, _ = step(call("run_analysis", design=True))
    assert r["ok"] and r["design"]["status"] in ("pass", "conditional", "fail", "incomplete")
    r, _, _ = step(call("read_results", document="design", path="/status"))
    assert r["ok"] and r["value"] == result({"tool_result": json.dumps(r)})["value"]
    for bad in ("garbage", json.dumps([1]), json.dumps({"tool": "fly"}), json.dumps({"tool": "inspect", "args": 3}),
                call("inspect", path=5), call("read_results", document="x")):
        r, _, _ = step(bad)
        assert r["ok"] is False and r["error"]
    r, terminated, info = step(call("submit", answer={"values": values}))
    # an agent that only uses the tools reads the reference values off the analysis: score 1
    assert terminated and rewards[-1] == 1.0 and all(x == 0.0 for x in rewards[:-1])
    assert info["score"]["score"] == 1.0 and info["reference"]["changeRecord"]["curveBall"] == CELLS["change"][1]
    assert info["steps_used"] == len(rewards)
    # after the end: no further reward, still terminated
    obs, reward, terminated, truncated, _ = env.step(call("inspect"))
    assert reward == 0.0 and terminated and "over" in result(obs)["error"]


def test_step_accepts_a_dict_action_too():
    env = make()
    env.reset(seed=0, options=options("change", task_type="predict"))
    obs, *_ = env.step({"tool": "inspect", "args": {"path": "/units"}})
    assert result(obs)["ok"]


def test_edit_is_a_validated_json_patch_on_the_working_copy():
    env = make(max_steps=60)
    obs, info = env.reset(seed=0, options=options("change", task_type="repair", mode="tool-allowed"))
    given = json.loads(obs["input_file"])
    first = next(iter(given["members"]))

    def do(action):
        return result(env.step(action)[0])

    # refused, nothing changed: a path that does not exist, an invalid model, malformed operations
    for ops in ([{"op": "replace", "path": "/members/nope/section", "value": "x"}],
                [{"op": "replace", "path": f"/members/{first}/section", "value": "no-such-section"}],
                [{"op": "remove", "path": "/nodes"}], [{"op": "move", "path": "/a"}], [{"op": "add", "path": "/x"}],
                [{"op": "add", "path": "no-slash", "value": 1}], [], "x", None,
                [{"op": "replace", "path": "/loadCases/" + next(iter(given["loadCases"])) + "/loads/99", "value": {}}]):
        r = do(call("edit", ops=ops))
        assert r["ok"] is False and r["applied"] == 0
        assert do(call("inspect"))["value"] == given
    r = do(call("edit", ops=[{"op": "replace", "path": f"/members/{first}/section", "value": "no-such-section"}]))
    assert r["problems"][0]["code"] == "ref.unknown-section"
    # several operations are applied together (an intermediate state need not be valid)
    sec = given["members"][first]["section"]
    r = do(call("edit", ops=[{"op": "add", "path": "/sections/copy", "value": given["sections"][sec]},
                             {"op": "replace", "path": f"/members/{first}/section", "value": "copy"},
                             {"op": "add", "path": f"/members/{first}/type", "value": "frame"}]))
    assert r["ok"] and r["applied"] == 3
    doc = do(call("inspect"))["value"]
    assert doc["members"][first]["section"] == "copy" and doc["sections"]["copy"] == given["sections"][sec]
    assert json.loads(env.step(call("inspect", path="/name"))[0]["input_file"]) == given     # the given file stays
    r = do(call("edit", ops=[{"op": "remove", "path": f"/members/{first}/type"}]))
    assert r["ok"] and "type" not in do(call("inspect", path=f"/members/{first}"))["value"]
    # array elements: add (insert / append), replace, remove
    lc = next(iter(given["loadCases"]))
    n = len(given["loadCases"][lc]["loads"])
    load = given["loadCases"][lc]["loads"][0]
    assert do(call("edit", ops=[{"op": "add", "path": f"/loadCases/{lc}/loads/-", "value": load}]))["ok"]
    assert len(do(call("inspect", path=f"/loadCases/{lc}/loads"))["value"]) == n + 1
    assert do(call("edit", ops=[{"op": "remove", "path": f"/loadCases/{lc}/loads/{n}"}]))["ok"]
    assert do(call("inspect", path=f"/loadCases/{lc}/loads"))["value"] == given["loadCases"][lc]["loads"]


def test_apply_edit_function():
    doc = {"a": {"b": [1, 2, 3]}, "c": 1}
    out = E.apply_edit(doc, [{"op": "add", "path": "/a/b/1", "value": 9}, {"op": "replace", "path": "/c", "value": 2},
                             {"op": "remove", "path": "/a/b/0"}, {"op": "add", "path": "/d~1e", "value": {}}])
    assert out == {"a": {"b": [9, 2, 3]}, "c": 2, "d/e": {}} and doc == {"a": {"b": [1, 2, 3]}, "c": 1}
    for ops in ([{"op": "remove", "path": "/zz"}], [{"op": "replace", "path": "/a/b/7", "value": 0}],
                [{"op": "add", "path": "/c/x", "value": 0}], [{"op": "add", "path": "/a/b/x", "value": 0}]):
        with pytest.raises(E.EditRefused):
            E.apply_edit(doc, ops)


def test_repair_in_the_environment_by_editing_the_working_copy():
    """The agent repairs the model with edit and submits the working copy; the environment re-solves it
    through the verification pipeline."""
    env = make(max_steps=10)
    obs, info = env.reset(seed=0, options=options("change", task_type="repair", mode="tool-allowed"))
    record = case("change")["record"]
    ops = []
    for c in reversed(record["changes"]):                        # undo the recorded change (the oracle's knowledge)
        ops.append({"op": "remove", "path": c["path"]} if c["op"] == "add" else
                   {"op": "add" if c["op"] == "remove" else "replace", "path": c["path"], "value": c["before"]})
    obs, reward, *_ = env.step(call("edit", ops=ops))
    assert result(obs)["ok"] and reward == 0.0
    obs, reward, terminated, truncated, info = env.step(call("submit"))
    assert terminated and reward == 1.0 and info["score"]["reasons"]["criteria"]["verified"]
    # submitting without a repair: nothing
    env.reset(seed=0, options=options("change", task_type="repair", mode="tool-allowed"))
    assert env.step(call("submit", answer={}))[1] == 0.0


def test_run_analysis_returns_a_refusal_as_an_observation():
    env = make()
    env.reset(seed=0, options=options("mechanism", task_type="diagnose", mode="tool-allowed"))
    obs, reward, terminated, truncated, info = env.step(call("run_analysis"))
    r = result(obs)
    assert r["ok"] is False and r["refusal"] == "unstable" and not terminated and reward == 0.0
    assert r["problems"][0]["code"] == "stability.mechanism" and r["diagnosis"]["nMechanisms"] >= 1
    assert result(env.step(call("read_results"))[0])["ok"] is False


def test_no_solver_mode_disables_the_analysis_tool():
    env = make()
    obs, info = env.reset(seed=0, options=options("change", task_type="predict", mode="no-solver"))
    assert info["mode"] == "no-solver" and "analysis tool is disabled" in obs["statement"]
    for action in (call("run_analysis"), call("run_analysis", design=True)):
        r = result(env.step(action)[0])
        assert r["ok"] is False and "disabled" in r["error"]
    r = result(env.step(call("read_results", path="/combinations"))[0])
    assert r["ok"] is False and "no analysis results" in r["error"]
    assert result(env.step(call("inspect", path="/units"))[0])["ok"]          # inspecting and editing still work
    assert result(env.step(call("edit", ops=[{"op": "replace", "path": "/name", "value": "x"}]))[0])["ok"]
    assert result(env.step(call("run_analysis"))[0])["ok"] is False


def test_truncation_on_the_step_budget():
    env = make(max_steps=3)
    obs, info = env.reset(seed=0, options=options("change", task_type="diagnose"))
    assert obs["steps_remaining"] == 3
    flags = []
    for k in range(3):
        obs, reward, terminated, truncated, info = env.step("x" if k == 1 else call("inspect", path="/units"))
        flags.append((reward, terminated, truncated, obs["steps_remaining"]))
    assert flags == [(0.0, False, False, 2), (0.0, False, False, 1), (0.0, False, True, 0)]
    assert "score" not in info and info["reference"]["changeRecord"]            # the episode is over: reference shown
    # a submit on the last allowed step terminates, it does not truncate
    env.reset(seed=0, options=options("change", task_type="diagnose"))
    env.step(call("inspect"))
    env.step(call("inspect"))
    obs, reward, terminated, truncated, info = env.step(call("submit", answer={"category": "supports"}))
    assert terminated and not truncated and reward == pytest.approx(S.DIAGNOSE_WEIGHTS["category"])


# -------------------------------------------------------------------------------------------- leak
def _forbidden(name):
    c = case(name)
    words = [c["key"]["curveBall"], c["record"]["target"], "changeRecord", "generatorSeed", "tiersPassed",
             "intendedModel", "whyItMatters", "baseHash"]
    return words


@pytest.mark.parametrize("name,task_type", [("change", "predict"), ("change", "diagnose"), ("change", "repair"),
                                            ("mechanism", "diagnose"), ("mechanism", "repair")])
def test_reference_does_not_leak_before_the_episode_ends(name, task_type):
    """No-solver mode (the analysis tool would legitimately print results): nothing of the hidden part --
    curve ball, change record, reference values, verification tiers, seeds -- is in an observation or in
    info before termination, whatever the agent calls."""
    env = make(max_steps=12)
    obs, info = env.reset(seed=3, options=options(name, task_type=task_type, mode="no-solver"))
    answers = env.hidden_answers()
    numbers = [repr(v["value"]) for v in answers["reference"].get("values", {}).values()]
    forbidden = _forbidden(name) + numbers
    if name == "mechanism":
        forbidden += ["nMechanisms", "\"kind\": \"diagnosis\"", "unrestrained"]
    actions = [call("inspect"), call("inspect", path="/members"), call("run_analysis"), call("read_results"),
               call("read_results", document="design"), "junk", call("edit", ops=[{"op": "replace", "path": "/name", "value": "y"}]),
               call("submit_answer"), call("inspect", path="/answers"), call("read_results", path="/hidden")]
    seen = [json.dumps([obs, info], default=str)]
    for action in actions:
        obs, reward, terminated, truncated, info = env.step(action)
        assert not terminated and not truncated and reward == 0.0
        seen.append(json.dumps([obs, info], default=str))
    for text in seen:
        for word in forbidden:
            assert word not in text, word
    # the info keys are the public ones
    assert set(info) <= {"task_id", "task_type", "mode", "family", "difficulty", "split", "steps_used", "max_steps",
                         "tools", "answer_format", "verification", "instances_skipped", "quantities", "criteria",
                         "catalogue"}
    assert info["verification"]["status"] == "verified" and "tiers" not in info["verification"]
    # what a task says about its verification is the same for a mechanism and for a stable model
    other = make_task(case("change" if name == "mechanism" else "mechanism"), task_type, "no-solver")[0] \
        if task_type != "predict" else None
    assert other is None or other["verification"] == info["verification"]
    # ... and when the episode ends the reference is given
    obs, reward, terminated, truncated, info = env.step(call("submit", answer={}))
    assert terminated and info["reference"] == answers and info["reference"]["verification"]["tiersPassed"]
    assert info["reference"]["verification"]["tiers"]["C"]["tolerances"]


def test_tool_allowed_info_does_not_leak_either():
    env = make()
    obs, info = env.reset(seed=3, options=options("change", task_type="predict", mode="tool-allowed"))
    answers = env.hidden_answers()
    text = json.dumps([obs, info])
    for word in _forbidden("change") + [repr(v["value"]) for v in answers["reference"]["values"].values()]:
        assert word not in text
    obs, _, _, _, info = env.step(call("run_analysis"))
    assert all(word not in json.dumps(info) for word in _forbidden("change"))


# ---------------------------------------------------------------------------------------- held out
def test_held_out_instances_are_refused_unless_the_environment_is_built_for_them():
    family, cb, seed, level = CELLS["held_out_seed"]
    assert splits.is_held_out(family, cb, seed)
    train, held = make(), make(split="held-out")
    with pytest.raises(ValueError, match="held-out"):
        train.reset(seed=0, options=options("held_out_seed"))
    obs, info = held.reset(seed=0, options=options("held_out_seed", task_type="diagnose"))
    assert info["split"] == "held-out"
    # a held-out COMBINATION is refused by the train environment whatever the seed
    hf, hc = splits.load_splits()["heldOutCombinations"][0]
    assert not splits.is_held_out(hf, "support_removed", 0) and splits.is_held_out(hf, hc, 0)
    for opts in ({"family": hf, "curve_ball": hc, "instance_seed": 0}, {"family": hf, "curve_ball": hc}):
        with pytest.raises(ValueError, match="held-out"):
            train.reset(seed=0, options=opts)
    # ... and the held-out environment refuses train instances
    with pytest.raises(ValueError, match="train"):
        held.reset(seed=0, options=options("change"))
    # random draws stay inside the split
    restrict = dict(families=("continuous_beam",), curve_balls=("support_fixed_to_pinned", "support_removed"),
                    difficulties=(2,), seed_range=5)
    for env, want in ((make(**restrict), False), (make(split="held-out", **restrict), True)):
        for s in range(3):
            env.reset(seed=s)
            key = env.hidden_answers()["key"]
            assert splits.is_held_out(key["family"], key["curveBall"], key["seed"]) is want


def test_an_instance_without_a_verified_reference_is_refused_with_the_log():
    env = make()
    with pytest.raises(TaskUnavailable) as e:
        env.reset(seed=0, options=options("design_only", task_type="diagnose"))
    assert "not-covered" in e.value.reason and e.value.log[0]["skips"]
    with pytest.raises(TaskUnavailable) as e:
        env.reset(seed=0, options=options("mechanism", task_type="predict"))
    assert "mechanism" in e.value.reason
    with pytest.raises(TaskUnavailable):
        env.reset(seed=0, options=options("change", task_type="design"))


# ------------------------------------------------------------------------------------------ agents
def test_oracle_agent_scores_one_and_baseline_does_not():
    env = gymnasium.make(E.ENV_ID, cache_dir=CACHE, max_steps=5)
    oracle, baseline, template = (load_agent(n) for n in ("oracle_agent", "random_agent", "llm_agent_template"))
    for name, types in (("change", ("predict", "diagnose", "repair")), ("mechanism", ("diagnose", "repair"))):
        for tt in types:
            for mode in ("tool-allowed", "no-solver"):
                opts = options(name, task_type=tt, mode=mode)
                ep = E.run_episode(env, oracle.make_agent(env), seed=0, options=opts)
                assert ep["reward"] == 1.0 and ep["terminated"] and ep["steps"] == 1
                ep = E.run_episode(env, baseline.make_agent(env), seed=0, options=opts)
                assert ep["reward"] < 0.5 and ep["terminated"]
    ep = E.run_episode(env, template.make_agent(env), seed=0, options=options("change", task_type="predict"))
    assert ep["reward"] == 0.0 and ep["terminated"]              # the stub submits an empty answer


def test_llm_template_takes_a_callable_and_survives_bad_replies():
    template = load_agent("llm_agent_template")
    env = make(max_steps=4)
    replies = iter(["I will look first.", "```json\n{\"tool\": \"inspect\", \"args\": {\"path\": \"/units\"}}\n```",
                    "{broken", "again nothing"])
    prompts = []

    def llm(prompt):
        prompts.append(prompt)
        return next(replies)

    ep = E.run_episode(env, template.make_agent(env, llm=llm), seed=0, options=options("change", task_type="diagnose"))
    assert ep["truncated"] and not ep["terminated"] and ep["reward"] == 0.0 and ep["steps"] == 4
    assert json.loads(ep["transcript"][1]["tool_result"])["ok"] is True
    assert "TASK (diagnose)" in prompts[0] and "TOOL RESULT" in prompts[1]
    assert template.extract_json('x {"a": {"b": "}"}} y') == {"a": {"b": "}"}} and template.extract_json("no") is None


def test_render_ansi():
    env = make(render_mode="ansi")
    assert env.render() is None
    env.reset(seed=0, options=options("change", task_type="predict"))
    assert env.render().startswith("TASK (predict)")
