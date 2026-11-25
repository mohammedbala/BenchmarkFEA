"""Baseline agent: looks at nothing that matters and submits a no-op answer.

    benchfea bench run-agent examples/agents/random_agent.py --dataset data --episodes 5 --seed 0

predict   every quantity 0
diagnose  a random category, a random member of the file, a random direction
repair    the model exactly as given (no repair)
design    the sections exactly as given

It shows the floor of the benchmark: what an agent scores without understanding the model.
"""
import json
import random


def make_agent(env, seed: int = 0):
    rng = random.Random(seed)

    def agent(observation, info):
        task_type = info["task_type"]
        if task_type == "predict":
            answer = {"values": {q["id"]: 0.0 for q in info["quantities"]}}
        elif task_type == "diagnose":
            model = json.loads(observation["input_file"])
            fmt = info["answer_format"]
            answer = {"category": rng.choice(fmt["category"]),
                      "entities": {"members": [rng.choice(sorted(model["members"]))]},
                      "direction": rng.choice(fmt["direction"]), "mechanism": None}
        elif task_type == "repair":
            answer = {}                       # the environment then submits the working copy: the model as given
        else:
            answer = {"sections": {}}
        return json.dumps({"tool": "submit", "args": {"answer": answer}})

    return agent
