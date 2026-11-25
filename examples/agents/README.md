# Example agents

Agents for the Gymnasium environment `BenchFEA-v0` (docs/BENCHMARK.md sections 5 and 6). Each file
defines `make_agent(env) -> agent`, where `agent(observation, info)` returns one tool call as JSON
text: `{"tool": "inspect" | "edit" | "run_analysis" | "read_results" | "submit", "args": {...}}`.

| file | what it is |
|---|---|
| `random_agent.py` | baseline: submits a no-op answer without looking at the model — the floor |
| `oracle_agent.py` | scripted: reads the hidden reference through a harness-side accessor and submits it — shows that every served task is solvable and scores 1.0; it is not an agent solving the task |
| `llm_agent_template.py` | template: `make_agent(env, llm=callable)` with `llm(prompt) -> text`; no API key, no network call; the default stub submits an empty answer |

```
benchfea bench run-agent examples/agents/oracle_agent.py --dataset data --episodes 5 --seed 0
benchfea bench run-agent examples/agents/random_agent.py --dataset data --episodes 5 --seed 0
benchfea bench run-agent examples/agents/llm_agent_template.py --dataset data --episodes 1
```
