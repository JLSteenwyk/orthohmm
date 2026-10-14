import json
import pytest
from types import SimpleNamespace

from benchmark_tools import benchmark_production


@pytest.mark.parametrize("cap", [None, 0, 400])
@pytest.mark.parametrize("support_mode", ["legacy", "bipartition"])
def test_main_passes_accuracy_profile_to_production_cli(tmp_path, monkeypatch, cap, support_mode):
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
        "--phylogeny", "reconcile",
        "--gene-tree-support-mode", support_mode,
        *([] if cap is None else ["--search-max-candidates-per-query", str(cap)]),
    ]) == 0

    command = captured["command"]
    profile_index = command.index("--accuracy_profile")
    assert command[profile_index + 1] == "high_sensitivity"
    if support_mode == "legacy":
        assert "--gene_tree_support_mode" not in command
    else:
        assert command[command.index("--gene_tree_support_mode") + 1] == support_mode
    if cap is None:
        assert "--search_max_candidates_per_query" not in command
    else:
        assert command[command.index("--search_max_candidates_per_query") + 1] == str(cap)
    payload = json.loads(result_json.read_text())
    assert payload["harness"]["command"] == command
