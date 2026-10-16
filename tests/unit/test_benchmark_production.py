import json
import pytest
from types import SimpleNamespace

from benchmark_tools import benchmark_production


@pytest.mark.parametrize("cap", [None, 0, 400])
@pytest.mark.parametrize("support_mode", ["legacy", "bipartition"])
@pytest.mark.parametrize("result_order", ["arrival", "canonical"])
@pytest.mark.parametrize("native_guard", [False, True])
@pytest.mark.parametrize("tie_policy", ["first_hit", "defer_ambiguous"])
def test_main_passes_accuracy_profile_to_production_cli(tmp_path, monkeypatch, cap, support_mode, result_order, native_guard, tie_policy):
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
        "--search-result-order", result_order,
        "--singleton-tie-policy", tie_policy,
        *(["--native-kernel-trace", str(tmp_path / "routes")] if native_guard else []),
        *([] if cap is None else ["--search-max-candidates-per-query", str(cap)]),
    ]) == 0

    command = captured["command"]
    if native_guard:
        assert command[1].endswith("benchmark_tools/strict_native_cli.py")
        assert command[2:4] == ["--trace-directory", str(tmp_path / "routes")]
    else:
        assert command[1:3] == ["-m", "orthohmm"]
    profile_index = command.index("--accuracy_profile")
    assert command[profile_index + 1] == "high_sensitivity"
    if tie_policy == "first_hit":
        assert "--singleton_tie_policy" not in command
    else:
        assert command[command.index("--singleton_tie_policy") + 1] == tie_policy
    if result_order=="arrival":
        assert "--search_result_order" not in command
    else:
        assert command[command.index("--search_result_order")+1]==result_order
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


@pytest.mark.parametrize("extra", [[], ["--accuracy-profile", "high_sensitivity"],
    ["--accuracy-profile", "high_sensitivity", "--search-result-order", "canonical", "--search-min-total-hits", "1"],
    ["--accuracy-profile", "high_sensitivity", "--search-result-order", "canonical", "--search-max-candidates-per-query", "0"]])
def test_union_benchmark_rejects_unsupported_settings(extra):
    with pytest.raises(SystemExit):
        benchmark_production.parse_args(["input", "output", "metrics.json", "--cpu", "2",
            "--search-candidate-policy", "protected_union", *extra])


def test_benchmark_forwards_opt_in_union_to_production(tmp_path, monkeypatch):
    source = tmp_path / "input"
    source.mkdir()
    (source / "s.faa").write_text(">gene\nACDE\n")
    captured = []
    def run(command, **kwargs):
        captured.extend(command)
        metrics = command[command.index("--metrics_json") + 1]
        with open(metrics, "w") as handle:
            json.dump({"schema_version": 1, "status": "complete"}, handle)
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(benchmark_production.subprocess, "run", run)
    monkeypatch.setattr(benchmark_production, "git_value", lambda *args: "")
    assert benchmark_production.main([str(source), str(tmp_path / "out"), str(tmp_path / "metrics.json"),
        "--cpu", "2", "--accuracy-profile", "high_sensitivity", "--search-result-order", "canonical",
        "--search-candidate-policy", "protected_union"]) == 0
    assert captured[captured.index("--search-candidate-policy") + 1] == "protected_union"
