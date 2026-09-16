import json

import pytest

from benchmark_tools.benchmark_production import file_record
from benchmark_tools.verify_simulation_histories import compare_histories, history_inventory, verify_generation


def native_fixture(root):
    for stage, parameter in (("T", "SpeciesTreeParameters.tsv"), ("G", "GenomeParameters.tsv")):
        (root / stage).mkdir(parents=True)
        (root / stage / "events.tsv").write_text("history\n")
        (root / stage / parameter).write_text("parameters\n")
    return root


def test_exact_history_comparison_excludes_only_parameter_files(tmp_path):
    a, b = [native_fixture(tmp_path / name) for name in ("a", "b")]
    (b / "G/GenomeParameters.tsv").write_text("different parameters")
    assert compare_histories(a, b)["matched"]
    (b / "G/events.tsv").write_text("different history")
    assert compare_histories(a, b)["differences"] == ["G/events.tsv"]
    (b / "T/extra.tsv").write_text("extra")
    assert compare_histories(a, b)["differences"] == ["G/events.tsv", "T/extra.tsv"]


def test_empty_stage_rejected(tmp_path):
    with pytest.raises(ValueError, match="No biological"):
        history_inventory(tmp_path)


def test_generation_completion_provenance_and_inventory(tmp_path):
    native = native_fixture(tmp_path / "native/run")
    evidence = tmp_path / "execution/run/status.json"
    evidence.parent.mkdir(parents=True)
    run = {"label": "run", "native_output": str(native), "commands": [{"stage": "T", "argv": ["command"]}]}
    status = {"label": "run", "status": "complete", "output_inventory_recorded": True,
              "provenance": {"manifest": {"sha256": "abc"}},
              "stages": [{"stage": "T", "command": ["command"], "status": "complete", "exit_code": 0}],
              "outputs": [file_record(p, tmp_path) for p in sorted(native.rglob("*")) if p.is_file()]}
    evidence.write_text(json.dumps(status))
    assert verify_generation(tmp_path, run, "abc") == native
    with pytest.raises(ValueError, match="different manifest"):
        verify_generation(tmp_path, run, "def")
    status["status"] = "running"
    evidence.write_text(json.dumps(status))
    with pytest.raises(ValueError, match="not verified complete"):
        verify_generation(tmp_path, run, "abc")
    status["status"] = "complete"
    status["outputs"] = []
    evidence.write_text(json.dumps(status))
    with pytest.raises(ValueError, match="absent"):
        verify_generation(tmp_path, run, "abc")
