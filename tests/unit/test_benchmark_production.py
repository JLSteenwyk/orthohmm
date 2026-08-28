import json
from types import SimpleNamespace

from benchmark_tools import benchmark_production


def test_main_passes_accuracy_profile_to_production_cli(tmp_path, monkeypatch):
    input_directory = tmp_path / "input"
    input_directory.mkdir()
    (input_directory / "species.faa").write_text(">gene\nACDE\n")
    output_directory = tmp_path / "output"
    result_json = tmp_path / "result.json"
    captured = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        metrics_path = command[command.index("--metrics_json") + 1]
        with open(metrics_path, "w") as handle:
            json.dump({"schema_version": 1, "status": "complete"}, handle)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(benchmark_production.subprocess, "run", fake_run)
    monkeypatch.setattr(benchmark_production, "git_value", lambda *args: "")

    assert benchmark_production.main([
        str(input_directory),
        str(output_directory),
        str(result_json),
        "--cpu", "2",
        "--accuracy-profile", "high_sensitivity",
    ]) == 0

    command = captured["command"]
    profile_index = command.index("--accuracy_profile")
    assert command[profile_index + 1] == "high_sensitivity"
    payload = json.loads(result_json.read_text())
    assert payload["harness"]["command"] == command
