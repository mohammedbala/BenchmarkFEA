"""Template: plug a language model into the environment through one callable.

    def my_llm(prompt: str) -> str:          # your code: call your model, return its text
        ...

    import gymnasium, benchfea.bench.env
    from llm_agent_template import make_agent
    env = gymnasium.make("BenchFEA-v0")
    agent = make_agent(env, llm=my_llm)
    from benchfea.bench.env import run_episode
    print(run_episode(env, agent, seed=1)["reward"])

``llm`` receives the whole conversation so far as one prompt and must return ONE tool call as JSON:
``{"tool": "inspect" | "edit" | "run_analysis" | "read_results" | "submit", "args": {...}}``.
This file contains no API key and makes no network call; without ``llm`` it uses a stub that
immediately submits an empty answer (score 0), so that the file runs as it is:

    benchfea bench run-agent examples/agents/llm_agent_template.py --dataset data --episodes 1
"""
import json

SYSTEM = """You are a structural engineer working with a frame analysis tool. Each turn, reply with exactly one JSON
object and nothing else: {"tool": NAME, "args": {...}}.
Tools:
  inspect       {"path": "/members/B1"}              read a value of the working copy of the model (JSON pointer)
  edit          {"ops": [{"op": "replace", "path": "/members/B1/section", "value": "..."}]}   change the working copy
  run_analysis  {"design": false}                    analyse the working copy (not available in mode no-solver)
  read_results  {"path": "/combinations/ULS1/displacements/N3", "document": "results"}
  submit        {"answer": {...}}                    ends the episode; the answer format is stated in the task
"""


def _stub_llm(prompt: str) -> str:
    return json.dumps({"tool": "submit", "args": {"answer": {}}})


def extract_json(text: str):
    """The first JSON object in a model reply (models often wrap it in prose or code fences)."""
    start = text.find("{")
    while start != -1:
        depth, in_string, escape = 0, False, False
        for k in range(start, len(text)):
            ch = text[k]
            if in_string:
                escape = (ch == "\\") and not escape
                if ch == '"' and not escape:
                    in_string = False
                continue
            if ch == '"':
                in_string = True
            elif ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    try:
                        return json.loads(text[start:k + 1])
                    except ValueError:
                        break
        start = text.find("{", start + 1)
    return None


def make_agent(env, llm=None, max_result_chars: int = 20000):
    llm = llm or _stub_llm
    history = []

    def agent(observation, info):
        if not history:                                       # first turn: the task and the model file
            history.append("TASK\n" + observation["statement"] + "\n\nMODEL FILE (JSON)\n" + observation["input_file"])
        else:
            history.append("TOOL RESULT\n" + observation["tool_result"][:max_result_chars])
        history.append(f"({observation['steps_remaining']} tool calls left)")
        reply = llm(SYSTEM + "\n\n" + "\n\n".join(history))
        history.append("YOUR CALL\n" + reply)
        call = extract_json(reply)
        # an unreadable reply is passed on as it is: the environment answers with an error observation
        return json.dumps(call) if call is not None else reply

    return agent
