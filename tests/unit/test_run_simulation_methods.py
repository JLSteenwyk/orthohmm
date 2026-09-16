import hashlib
import json
import os
from pathlib import Path
import sys

import pytest

from benchmark_tools.benchmark_production import file_record
from benchmark_tools.run_simulation_methods import copy_inputs, execute, read_frozen, verify_inputs


def test_frozen_hash_gate(tmp_path):
    path = tmp_path / "manifest.json"
    path.write_text('{"ok": true}')
    assert read_frozen(path, hashlib.sha256(path.read_bytes()).hexdigest()) == {"ok": True}
    with pytest.raises(ValueError, match="changed"):
        read_frozen(path, "0" * 64)


def test_copy_checks_content_and_refuses_existing_target(tmp_path):
    source = tmp_path / "input"
    source.mkdir()
    path = source / "a.fasta"
    path.write_text(">a\nAAA\n")
    record = dict(file_record(path, source), absolute_path=str(path))
    target = tmp_path / "copy"
    copy_inputs(source, target, [record])
    assert (target / path.name).read_bytes() == path.read_bytes()
    with pytest.raises(FileExistsError):
        copy_inputs(source, target, [record])
    path.write_text("changed")
    with pytest.raises(ValueError, match="changed"):
        copy_inputs(source, tmp_path / "copy2", [record])


def test_failed_method_preserved_and_next_method_executed(tmp_path):
    output = tmp_path / "success"
    dataset = {"label": "fixture", "methods": {
        "first": {"output": str(tmp_path / "failure"), "argv": [sys.executable, "-c", "raise SystemExit(7)"]},
        "second": {"output": str(output), "argv": [sys.executable, "-c",
            "from pathlib import Path; import sys; p=Path(sys.argv[1]); p.mkdir(); (p/'out').write_text('result')", str(output)]}}}
    evidence = tmp_path / "evidence"
    result = execute(dataset, ["first", "second"], os.environ.copy(), evidence, {"status": "ready", "inputs": []}, {})
    assert result["failed_methods"] == ["first"]
    assert result["methods"]["first"]["exit_code"] == 7
    assert result["methods"]["second"]["status"] == "process_succeeded"
    assert result["methods"]["second"]["outputs"][0]["sha256"] == hashlib.sha256(b"result").hexdigest()
    assert result["status"] == "finished_pending_native_validation"
    assert not result["native_outputs_validated"] and not result["accuracy_evaluated"]
    assert json.loads((evidence / "status.json").read_text()) == result
    with pytest.raises(FileExistsError):
        execute(dataset, ["first", "second"], os.environ.copy(), evidence, {"status": "ready"}, {})


def test_inapplicable_runs_no_command(tmp_path):
    dataset = {"label": "fixture", "methods": {"first": {"output": str(tmp_path / "out"), "argv": ["not-an-executable"]}}}
    result = execute(dataset, ["first"], {}, tmp_path / "evidence", {"status": "inapplicable", "reason": "fixture"}, {})
    assert result["status"] == "inapplicable" and result["methods"] == {}


def test_missing_generated_species_file_rejected(tmp_path, monkeypatch):
    inputs = tmp_path / "prepared/input"
    inputs.mkdir(parents=True)
    truth = inputs.parent / "truth.json"
    truth.write_text("{}")
    for name in ("a", "b", "c"):
        (inputs / f"{name}.fasta").write_text(f">{name}\nAAA\n")
    records = [file_record(p, tmp_path) for p in [truth, *inputs.iterdir()]]
    evidence = tmp_path / "execution/base/status.json"
    evidence.parent.mkdir(parents=True)
    evidence.write_text(json.dumps({"outputs": records}))
    dataset = {"parent": "base", "input": str(inputs), "truth": str(truth)}
    generation = {"simulation_runs": [{"label": "base"}, {"label": "divergent"}],
                  "history_equivalence_checks": [{"first": "base", "second": "divergent"}]}
    monkeypatch.setattr("benchmark_tools.run_simulation_methods.verify_generation", lambda *a: tmp_path)
    monkeypatch.setattr("benchmark_tools.run_simulation_methods.compare_histories", lambda *a: {"matched": True})
    assert verify_inputs(dataset, generation, tmp_path)["status"] == "ready"
    (inputs / "c.fasta").unlink()
    with pytest.raises(ValueError, match="file set changed"):
        verify_inputs(dataset, generation, tmp_path)
