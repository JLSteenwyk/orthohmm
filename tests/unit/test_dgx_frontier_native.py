import copy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import run_dgx_frontier_native_smoke as launch
from benchmark_tools.measure_native_frontier_step import measure


ROOT = Path(__file__).resolve().parents[2] / "benchmark_tools"


def test_transferred_recipe_matches_frozen_sources():
    import hashlib

    manifest = json.loads((ROOT / "results/dgx_frontier_native_recipe_v1_20260918.json").read_text())
    assert manifest["external_symlinks"] == []
    files = [row for row in manifest["records"] if row["kind"] == "file"]
    assert len(files) == 33
    for row in files:
        path = ROOT / Path(row["path"]).name
        frozen = {"probe_cgroup_frontier.py": "probe_cgroup_frontier_4ce76668.py.txt",
                  "measure_native_frontier_step.py": "measure_native_frontier_step_9e89534.py.txt"}
        if path.name in frozen:
            path = Path(__file__).parent / "fixtures" / frozen[path.name]
        if not path.exists():
            path = ROOT / "results" / path.name
        data = path.read_bytes()
        assert len(data) == row["bytes"]
        assert hashlib.sha256(data).hexdigest() == row["sha256"]


def test_retained_native_results_replay_without_relaxing_flags():
    from benchmark_tools.measure_native_frontier_step import evaluate

    report = json.loads((ROOT / "results/dgx_frontier_native_smokes_21831.json").read_text())
    assert report["scientific_timings_admitted"] == 0
    assert report["publication_ready"] is False
    flags = []
    for row in report["runs"]:
        measured = row["verification"]["measurement"]
        assert evaluate(measured["points"], measured["native"], measured["job_id"]) == measured["screening"]
        assert measured["native"]["exit_code"] == 0
        assert row["output_validation"]["input_genes"] == 645
        screen = measured["screening"]["original_threshold_screen"]
        flags.append(screen["flagged_intervals"])
        assert screen["whole_command_screen"]["screen_passed"] is True
        for token in ("JobState=COMPLETED", "Restarts=0", "CPUs/Task=20", "MinMemoryNode=96G", "OverSubscribe=NO"):
            assert token in row["scheduler"].split()
    assert flags == [[1], [1, 3], []]
    satellite = report["runs"][1]["verification"]["measurement"]["screening"]
    assert satellite["frontier_intervals"][3]["outside_target_frontier_cpu_s"] == pytest.approx(.001811)
    assert satellite["frontier_intervals"][3]["root_minus_frontier_cpu_s"] == pytest.approx(.377656)
    assert any(row["root_minus_frontier_cpu_s"] < 0 for row in satellite["frontier_intervals"])


@pytest.mark.parametrize("old,new", [
    ("run_dgx_hierarchy_quiet_smoke.py", "run_dgx_frontier_native_smoke.py"),
    ("run_dgx_hierarchy_quiet_smoke.sh", "run_dgx_frontier_native_smoke.sh"),
    ("audit_hierarchy_quiet_smokes.py", "audit_frontier_native_smokes.py"),
])
def test_only_collector_paths_and_reporting_differ(old, new):
    expected = (ROOT / old).read_text()
    replacements = {
        "hierarchy_quiet_smoke_v1": "frontier_native_smoke_v1",
        "hierarchy_quiet_recipe_v1": "frontier_native_recipe_v1",
        "run_dgx_hierarchy_quiet_smoke": "run_dgx_frontier_native_smoke",
        "three_hierarchy_quiet_smokes_validated": "three_frontier_native_smokes_validated",
        "from measure_native_hierarchy_step import measure": "from measure_native_frontier_step import measure",
        "from benchmark_tools.measure_native_hierarchy_step import evaluate, interval_point":
            "from benchmark_tools.measure_native_frontier_step import evaluate, interval_point",
        "measurement/hierarchy_report.json": "measurement/frontier_report.json",
        "complete-command hierarchy observation": "complete-command frontier observation",
        "Replay complete-command hierarchy smokes": "Replay complete-command frontier smokes",
    }
    for before, after in replacements.items():
        expected = expected.replace(before, after)
    assert (ROOT / new).read_text() == expected


def test_relocation_is_boundary_aware_and_does_not_mutate():
    prefix = str(launch.ROOT / "launcher_smoke_v1")
    original = {"values": [prefix + "/run_00", prefix + "_other/run_00", 3]}
    saved = copy.deepcopy(original)
    result = launch.relocate(original)
    assert result["values"][0] == str(launch.ROOT / "frontier_native_smoke_v1/run_00")
    assert result["values"][1:] == original["values"][1:]
    assert original == saved


@pytest.fixture
def environment(monkeypatch):
    monkeypatch.setattr(launch.os, "uname", lambda: SimpleNamespace(nodename="spark-7ff0"))
    monkeypatch.setattr(launch.sys, "dont_write_bytecode", True)
    for key, value in {"SLURM_CPUS_PER_TASK": "20", "SLURM_MEM_PER_NODE": "98304",
                       "SLURM_JOB_ID": "123", "PYTHONHASHSEED": "0"}.items():
        monkeypatch.setenv(key, value)
    for key in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setattr(Path, "exists", lambda self: False)
    monkeypatch.setattr(launch.os, "chdir", lambda path: None)


@pytest.mark.parametrize("index", [0, 1, 2])
def test_launch_preserves_frozen_native_command_and_pins_recipe(environment, monkeypatch, index):
    spec_path = ROOT / "results/dgx_launcher_smoke_spec_20260917.json"
    spec = json.loads(spec_path.read_text())
    calls = []
    def measured(*args, **kwargs):
        calls.append((args, kwargs))
        return {"status": "command_exited_zero"}
    monkeypatch.setattr(launch, "measure_run", measured)
    launch.launch(spec_path, index, Path("/frozen/recipe.json"), "a" * 64)
    args, kwargs = calls[0]
    assert args[0] == launch.relocate(spec["runs"][index])
    assert args[1] == spec["orders"][0]
    assert args[2] == [(row["path"], row["sha256"]) for row in spec["runtime_manifests"]] + [("/frozen/recipe.json", "a" * 64)]
    # The script imports by its local module name; compare source paths, not module identities.
    assert Path(args[4].__code__.co_filename).resolve() == Path(measure.__code__.co_filename).resolve()
    assert args[5] == 123
    assert kwargs == {"timeout_s": 900}


@pytest.mark.parametrize("key,value", [("SLURM_CPUS_PER_TASK", "19"),
                                      ("SLURM_MEM_PER_NODE", "65536"),
                                      ("PYTHONHASHSEED", "1"),
                                      ("LD_PRELOAD", "/unexpected.so")])
def test_incompatible_environment_rejected_before_execution(environment, monkeypatch, key, value):
    monkeypatch.setenv(key, value)
    with pytest.raises(ValueError):
        launch.launch(Path("/not/read"), 0, Path("/recipe"), "a" * 64)


@pytest.mark.parametrize("index", [-1, 3])
def test_out_of_range_method_rejected(index):
    with pytest.raises(ValueError, match="three frozen"):
        launch.launch(Path("/not/read"), index, Path("/recipe"), "a" * 64)
