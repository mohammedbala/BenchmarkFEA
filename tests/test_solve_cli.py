"""``benchfea run``: writes the results document, prints a summary, refuses with a diagnosis."""
import json

import pytest

from benchfea import io
from benchfea.cli import main

from solve_reference import make


def write(tmp_path, doc, name="model.json"):
    path = tmp_path / name
    io.dump(doc, path)
    return path


def cantilever(**kw):
    return make({"a": (0, 0, 0), "b": (3000, 0, 0)}, {"m": ("a", "b")}, {"a": "fixed"},
                {"c": [{"type": "nodal", "node": "b", "Fz": -1000.0}]}, **kw)


def test_run_writes_the_results_document_and_prints_a_summary(tmp_path, capsys, fixtures):
    out = tmp_path / "beam.results.json"
    assert main(["run", str(fixtures / "valid" / "beam.json"), "-o", str(out)]) == 0
    text = capsys.readouterr().out
    assert "1 combination(s); units kN, m" in text and "invariants passed" in text and str(out) in text
    assert "sum of reactions (0, 101.25, 0) kN" in text                      # 1.35 * 12.5 kN/m * 6 m
    doc = json.loads(out.read_text())
    assert doc["format"] == "benchfea-results" and doc["formatVersion"] == "1.0.0" and doc["units"] == {"force": "kN", "length": "m"}
    c = doc["combinations"]["ULS"]
    assert c["reactions"]["a"]["Fy"] == pytest.approx(50.625) and c["invariants"]["passed"] is True
    assert c["members"]["m1"]["stations"]["Mz"][5] == pytest.approx(1.35 * 12.5 * 36 / 8)      # statics: w L^2 / 8
    assert doc["defaultsUsed"] and isinstance(doc["warnings"], list) and "ULS" in doc["envelopes"]["m1"]["maxCombination"]["Mz"]


def test_run_default_output_name_yaml_input_and_json_summary(tmp_path, capsys, fixtures):
    src = tmp_path / "beam.yaml"
    src.write_text((fixtures / "valid" / "beam.yaml").read_text())
    assert main(["run", str(src), "--json"]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["ok"] is True and summary["output"] == str(tmp_path / "beam.results.json")
    assert summary["combinations"]["ULS"]["invariantsPassed"] is True
    assert summary["combinations"]["ULS"]["maxStationDisplacement"]["member"] == "m1"
    assert (tmp_path / "beam.results.json").exists()
    # same numbers as the JSON twin of the file
    out2 = tmp_path / "twin.json"
    assert main(["run", str(fixtures / "valid" / "beam.json"), "-o", str(out2)]) == 0
    a, b = json.loads((tmp_path / "beam.results.json").read_text()), json.loads(out2.read_text())
    assert a["combinations"] == b["combinations"]


def test_run_refuses_a_mechanism_with_the_diagnosis_and_writes_nothing(tmp_path, capsys):
    doc = cantilever()
    doc["supports"] = {"a": {"type": "pinned"}}
    path = write(tmp_path, doc)
    assert main(["run", str(path)]) == 1
    err = capsys.readouterr().err
    assert "stability.mechanism" in err and "UNSTABLE" in err and "rigid-body rotation" in err
    assert not (tmp_path / "model.results.json").exists()
    assert main(["run", str(path), "--json"]) == 1
    payload = json.loads(capsys.readouterr().out)
    assert payload["ok"] is False and payload["refusal"] == "unstable"
    assert payload["diagnosis"]["nMechanisms"] == 3 and payload["diagnosis"]["modes"][0]["kind"] == "rigid-body"


def test_run_refuses_invalid_and_unsupported_models(tmp_path, capsys, fixtures, monkeypatch):
    doc = cantilever()
    doc["members"]["m"]["section"] = "nope"
    assert main(["run", str(write(tmp_path, doc))]) == 1
    assert "ref.unknown-section" in capsys.readouterr().err
    # second-order analysis is implemented now; an analysis type without a registered solver is refused
    from benchfea.solve import ANALYSIS_TYPES
    monkeypatch.delitem(ANALYSIS_TYPES, "secondOrder")
    assert main(["run", str(fixtures / "valid" / "frame_kip_in.json"), "-o", str(tmp_path / "x.json")]) == 1
    err = capsys.readouterr().err
    assert "analysis.not-implemented" in err and "REFUSED" in err and not (tmp_path / "x.json").exists()


def test_run_flagged_results_need_an_explicit_switch(tmp_path, capsys, monkeypatch):
    from benchfea.solve.element import MemberElement
    original = MemberElement._fixed_end
    monkeypatch.setattr(MemberElement, "_fixed_end", lambda self, loads: tuple(1.01 * x for x in original(self, loads)))
    doc = cantilever()
    doc["loadCases"]["c"]["loads"].append({"type": "memberDistributed", "member": "m", "direction": "Z", "w1": -1.0})
    path = write(tmp_path, doc)
    assert main(["run", str(path)]) == 1
    assert "solve.invariant-failed" in capsys.readouterr().err and not (tmp_path / "model.results.json").exists()
    assert main(["run", str(path), "--allow-flagged"]) == 3
    assert "FLAGGED" in capsys.readouterr().out
    doc_out = json.loads((tmp_path / "model.results.json").read_text())
    assert doc_out["combinations"]["c"]["invariants"]["passed"] is False
    assert any("FLAGGED" in w for w in doc_out["warnings"])
