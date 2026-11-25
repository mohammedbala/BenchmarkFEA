"""Scripted "oracle" agent: reads the hidden reference of the task and submits it.

    benchfea bench run-agent examples/agents/oracle_agent.py --dataset data --episodes 5 --seed 0

This is NOT an agent solving the task: it uses ``env.unwrapped.hidden_answers()``, a harness-side
accessor that no observation or info gives.  It exists to show that every task the environment serves
is solvable and that the perfect answer scores 1.0 through the same scoring path as any other answer
(repairs and designs are re-solved through the verification pipeline).
"""
import json

from benchfea.bench.scoring import reference_answer


def make_agent(env):
    base = env.unwrapped

    def agent(observation, info):
        answers = base.hidden_answers()
        task = {"taskType": info["task_type"]}
        return json.dumps({"tool": "submit", "args": {"answer": reference_answer(task, answers)}})

    return agent
