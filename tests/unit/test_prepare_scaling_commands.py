import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import benchmark_production as harness
from benchmark_tools.prepare_scaling_commands import configurations, native_orthohmm, prepare
from benchmark_tools.prepare_scaling_inputs import planned_runs
from benchmark_tools.prepare_simulation_methods import commands


@pytest.mark.parametrize("method", ["orthohmm_high_sensitivity", "orthohmm_satellite_v2"])
@pytest.mark.parametrize("symlinked_inputs", [False, True])
def test_native_command_matches_actual_harness_translation(tmp_path, monkeypatch, method, symlinked_inputs):
    inputs = tmp_path / "inputs"
    if symlinked_inputs:
        physical = tmp_path / "physical_inputs"
        physical.mkdir()
        inputs.symlink_to(physical, target_is_directory=True)
    else:
        inputs.mkdir()
    core = tmp_path / "core"
    config = commands({"input": str(inputs)}, tmp_path / "outputs", core,
                      sys.executable, tmp_path / "orthofinder")[method]
    config["argv"][config["argv"].index("--cpu") + 1] = "32"
    captured = []
    monkeypatch.setattr(harness, "git_value", lambda *a: "fixture")
    monkeypatch.setattr(harness.subprocess, "run", lambda argv, **kw: (captured.append(argv) or SimpleNamespace(returncode=0)))
    assert harness.main(config["argv"][2:]) == 0
    assert captured == [native_orthohmm(config)]
    assert "benchmark_production.py" not in " ".join(captured[0])
    assert captured[0][captured[0].index("--threads_per_worker") + 1] == "4"


def fixtures(tmp_path):
    inputs = {"planned_runs": planned_runs(), "datasets": [
        {"proteomes": n, "input_directory": str(tmp_path / str(n))} for n in (4, 8, 12)]}
    baseline = {"core_root": str(tmp_path / "core"), "tool_entrypoints": {
        "orthohmm_python": {"absolute_path": sys.executable},
        "orthofinder": {"absolute_path": str(tmp_path / "orthofinder")}}}
    return inputs, baseline


def test_all_runs_fresh_balanced_and_resource_only_changes(tmp_path):
    inputs, baseline = fixtures(tmp_path)
    runs = configurations(inputs, baseline, tmp_path / "runs")
    assert len(runs) == 27
    assert [{k: r[k] for k in ("index", "repeat", "proteomes", "method")} for r in runs] == planned_runs()
    assert len({r["configuration"]["output"] for r in runs}) == 27
    assert len({r["measurement_directory"] for r in runs}) == 27
    assert not (tmp_path / "runs").exists()
    for run in runs:
        argv = run["native_argv"]
        if run["native_method"] == "orthofinder_full":
            assert argv[argv.index("-t") + 1] == argv[argv.index("-a") + 1] == "32"
            assert "-s" not in argv and "-ft" not in argv and "-b" not in argv
            assert run["configuration"]["copy_inputs_to"].startswith(run["configuration"]["output"])
        else:
            assert argv[argv.index("-c") + 1] == "32"
            assert argv[argv.index("--threads_per_worker") + 1] == "4"
            assert argv[argv.index("--accuracy_profile") + 1] == "high_sensitivity"
            if run["native_method"].endswith("satellite_v2"):
                assert argv[argv.index("--species_tree_mode") + 1] == "infer"
                assert argv[argv.index("--phylogeny_candidates") + 1] == "satellite_v2"


@pytest.mark.parametrize("change", ["missing_run", "reorder", "missing_size"])
def test_changed_panel_rejected(tmp_path, change):
    inputs, baseline = fixtures(tmp_path)
    if change == "missing_run":
        inputs["planned_runs"].pop()
    elif change == "reorder":
        inputs["planned_runs"].reverse()
    else:
        inputs["datasets"].pop()
    with pytest.raises(ValueError):
        configurations(inputs, baseline, tmp_path / "runs")


def test_existing_destination_rejected(tmp_path):
    with pytest.raises(FileExistsError):
        prepare(tmp_path, tmp_path, tmp_path / "manifest.json")


def test_dgx_resource_change_preserves_all_other_settings(tmp_path):
    inputs, baseline = fixtures(tmp_path)
    original = configurations(inputs, baseline, tmp_path / "runs")
    dgx = configurations(inputs, baseline, tmp_path / "runs", cpu_count=20)
    assert len(dgx) == len(original) == 27
    for before, after in zip(original, dgx):
        # Normalize only declared allocation flags, then compare the entire record.
        changes = []
        for key, flags in (("native_argv", ("-t", "-a") if before["native_method"] == "orthofinder_full" else ("-c",)),
                           ("configuration", ("-t", "-a") if before["native_method"] == "orthofinder_full" else ("--cpu",))):
            argv = after[key]["argv"] if key == "configuration" else after[key]
            for flag in flags:
                assert argv[argv.index(flag) + 1] == "20"
                changes.append((argv, argv.index(flag) + 1))
        for argv, index in changes:
            argv[index] = "32"
        assert after == before


@pytest.mark.parametrize("cpus", [0, -1, True, 20.0, "20", None])
def test_invalid_cpu_count_rejected(tmp_path, cpus):
    inputs, baseline = fixtures(tmp_path)
    with pytest.raises(ValueError, match="CPU allocation"):
        configurations(inputs, baseline, tmp_path / "runs", cpu_count=cpus)


@pytest.mark.parametrize("method, supplied_tree", [
    ("orthohmm_high_sensitivity", False),
    ("orthohmm_satellite_v2", False),
    ("orthohmm_satellite_v2", True),
])
def test_remote_native_paths_preserve_lexical_target(tmp_path, method, supplied_tree):
    physical = tmp_path / "physical"
    physical.mkdir()
    remote = tmp_path / "remote"
    remote.symlink_to(physical, target_is_directory=True)
    config = commands({"input": str(remote / "inputs")}, remote / "outputs", remote / "core",
                      remote / "python", remote / "orthofinder")[method]
    if supplied_tree:
        config["argv"] += ["--species-tree", str(remote / "species.nwk")]
        config["argv"][config["argv"].index("--species-tree-mode") + 1] = "supplied"
    local = native_orthohmm(config)
    target = native_orthohmm(config, resolve_paths=False)
    assert local[0] == target[0] == str(remote / "python")
    indices = [3, target.index("-o") + 1, target.index("--metrics_json") + 1]
    if supplied_tree:
        indices.append(target.index("--species_tree") + 1)
    for index in indices:
        assert Path(target[index]).is_relative_to(remote)
        assert local[index] == str(physical / Path(target[index]).relative_to(remote))
        local[index] = target[index]
    assert local == target
    assert not (physical / "outputs").exists()


@pytest.mark.parametrize("field", ["python", "input", "output", "metrics", "species_tree"])
@pytest.mark.parametrize("path", ["relative/path", "/remote/../path"])
def test_remote_native_paths_reject_ambiguity(tmp_path, field, path):
    config = commands({"input": str(tmp_path / "inputs")}, tmp_path / "outputs", tmp_path / "core",
                      tmp_path / "python", tmp_path / "orthofinder")["orthohmm_satellite_v2"]
    if field == "species_tree":
        config["argv"] += ["--species-tree", path]
    else:
        config["argv"][{"python": 0, "input": 2, "output": 3, "metrics": 4}[field]] = path
    with pytest.raises(ValueError, match="absolute without parent traversal"):
        native_orthohmm(config, resolve_paths=False)


@pytest.mark.parametrize("mode", [None, 0, 1, "False"])
def test_invalid_path_resolution_mode_rejected(tmp_path, mode):
    inputs, baseline = fixtures(tmp_path)
    config = commands({"input": str(tmp_path / "inputs")}, tmp_path / "outputs", tmp_path / "core",
                      tmp_path / "python", tmp_path / "orthofinder")["orthohmm_high_sensitivity"]
    with pytest.raises(ValueError, match="boolean path-resolution"):
        native_orthohmm(config, resolve_paths=mode)
    with pytest.raises(ValueError, match="boolean path-resolution"):
        configurations(inputs, baseline, tmp_path / "runs", resolve_orthohmm_paths=mode)
