import copy
import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import qfo_private_phylogeny_environment as environment
from benchmark_tools import run_qfo_private_phylogeny_control as driver
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.prepare_orthobench_factorial import plan_cells


def deployment_fixture():
    legacy = {"environments": {"orthofinder": {"packages": {"orthohmm": "0.5.0", "numpy": "x"}},
        "orthohmm": {"packages": {"numpy": "2.2.6", "unrelated": "1"}}},
        "tool_entrypoints": {"orthohmm_python": {"absolute_path": "old"}, "mafft": "frozen"},
        "core_sources": ["immutable"], "environment_overrides": {"OMP_NUM_THREADS": "1"}}
    ancestral = copy.deepcopy(legacy)
    ancestral["environments"]["orthofinder"]["packages"].pop("orthohmm")
    ancestral["threadripper_amendment"] = {"scope": "cwd correction"}
    python = dict(path="private", bytes=1, sha256="pin")
    deployment = {"new_interpreter": python}
    amended = copy.deepcopy(ancestral)
    amended["tool_entrypoints"]["orthohmm_python"] = dict(path="python", absolute_path="private", bytes=1, sha256="pin")
    amended["environments"]["orthohmm"] = {"packages": {"numpy": "2.2.6"}}
    return legacy, ancestral, amended, deployment


@pytest.mark.parametrize("problem", [None, "metadata", "ancestor_source", "ancestor_tool", "source", "tool", "threads", "orthofinder", "interpreter"])
def test_exact_deployment_difference(problem):
    legacy, ancestral, amended, deployment = deployment_fixture()
    before = copy.deepcopy((legacy, ancestral, amended, deployment))
    if problem == "metadata":
        legacy["environments"]["orthofinder"]["packages"]["orthohmm"] = "changed"
    if problem == "ancestor_source":
        ancestral["core_sources"] = []
    if problem == "ancestor_tool":
        ancestral["tool_entrypoints"]["mafft"] = "changed"
    if problem == "source":
        amended["core_sources"] = []
    if problem == "tool":
        amended["tool_entrypoints"]["mafft"] = "changed"
    if problem == "threads":
        amended["environment_overrides"]["OMP_NUM_THREADS"] = "32"
    if problem == "orthofinder":
        amended["environments"]["orthofinder"]["packages"]["numpy"] = "changed"
    if problem == "interpreter":
        amended["tool_entrypoints"]["orthohmm_python"]["absolute_path"] = "shared"
    if problem:
        with pytest.raises(ValueError):
            environment.deployment_difference(legacy, ancestral, amended, deployment)
    else:
        environment.deployment_difference(legacy, ancestral, amended, deployment)
        assert (legacy, ancestral, amended, deployment) == before


def original_cell(root):
    return next(cell for cell in plan_cells(root / "baseline", root / "fastas",
        root / "prepared/benchmark_tools/replay_phylogeny.py", 32) if cell["label"] == "p1_c1_r1")


@pytest.mark.parametrize("problem", [None, "label", "supplied", "tree", "checkpoint", "scoring", "cpu", "duplicate"])
def test_control_preserves_science(tmp_path, monkeypatch, problem):
    from benchmark_tools import run_qfo_factorial_cell as commands

    cell = original_cell(tmp_path)
    argv = cell["argv"]
    if problem == "label":
        cell["label"] = "p0_c1_r1"
    if problem == "supplied":
        argv[argv.index("--species-tree-mode") + 1] = "supplied"
    if problem == "tree":
        argv.extend(["--species-tree", "reference-tree"])
    if problem == "checkpoint":
        argv.extend(["--checkpoint-source", "other"])
    if problem == "scoring":
        argv.extend(["--official-benchmark", "labels"])
    if problem == "cpu":
        argv[argv.index("--cpu") + 1] = "64"
    if problem == "duplicate":
        argv.extend(["--json", "duplicate"])
    before = copy.deepcopy(cell)
    monkeypatch.setattr(commands, "native_command", lambda row, *args: (list(row["argv"]), []))
    private = {"tool_entrypoints": {"orthohmm_python": {"absolute_path": "private-python"}}}
    if problem:
        with pytest.raises(ValueError):
            environment.control_command(cell, tmp_path / "launcher", tmp_path / "prepared", private, tmp_path / "control")
    else:
        actual, _ = environment.control_command(cell, tmp_path / "launcher", tmp_path / "prepared", private, tmp_path / "control")
        expected = list(before["argv"])
        checkpoint = expected[expected.index("--output-directory") + 1]
        expected[0] = "private-python"
        expected[expected.index("--output-directory") + 1] = str(tmp_path / "control/output")
        expected[expected.index("--json") + 1] = str(tmp_path / "control/metrics.json")
        expected.extend(["--checkpoint-source", checkpoint])
        assert actual == expected and cell == before


@pytest.mark.parametrize("changed", [None, *driver.NATIVE_FILES])
def test_complete_native_comparison(tmp_path, changed):
    old, new = tmp_path / "old", tmp_path / "new"
    old.mkdir()
    new.mkdir()
    for name in driver.NATIVE_FILES:
        (old / name).write_text(name)
        (new / name).write_text(name + ("changed" if changed == name else ""))
    result = driver.compare_outputs(old, new)
    assert set(result) == set(driver.NATIVE_FILES)
    assert [name for name, item in result.items() if not item["byte_equal"]] == ([] if changed is None else [changed])


def allocation(monkeypatch):
    for key, value in dict(SLURM_JOB_ID="99999", SLURM_CPUS_PER_TASK="32", SLURM_MEM_PER_NODE="196608",
                          SLURM_JOB_NODELIST="bizon").items():
        monkeypatch.setenv(key, value)
    monkeypatch.delenv("SLURM_ARRAY_TASK_ID", raising=False)
    monkeypatch.setattr(driver.sys, "executable", "/home/bizon/anaconda3/bin/python")


@pytest.mark.parametrize("key,value", [("SLURM_JOB_ID", ""), ("SLURM_CPUS_PER_TASK", "2"),
    ("SLURM_MEM_PER_NODE", "65536"), ("SLURM_JOB_NODELIST", "dgx"), ("SLURM_ARRAY_TASK_ID", "1")])
def test_allocation_before_reads(tmp_path, monkeypatch, key, value):
    allocation(monkeypatch)
    monkeypatch.setenv(key, value)
    with pytest.raises(ValueError, match="standalone"):
        driver.run(tmp_path, "unused")


@pytest.mark.parametrize("outcome", ["success", "parity", "failed", "exception", "changed", "lookup", "protocol"])
def test_control_execution_is_separate_and_fail_closed(tmp_path, monkeypatch, outcome):
    from benchmark_tools import run_simulation_methods as execution_module

    allocation(monkeypatch)
    launcher = tmp_path / "launcher"
    launcher.mkdir()
    cell = original_cell(tmp_path)
    baseline = Path(cell["argv"][cell["argv"].index("--output-directory") + 1]) / "orthohmm_phylogeny"
    baseline.mkdir(parents=True)
    for name in driver.NATIVE_FILES:
        (baseline / name).write_text(name)
    before = copy.deepcopy(cell)
    protocol = tmp_path / driver.PROTOCOL
    protocol.parent.mkdir(parents=True)
    protocol.write_text("protocol")
    verified = {"launcher": str(launcher), "prepared": str(tmp_path / "prepared"), "original": cell,
        "manifest": {"input_fastas": [], "environment_overrides": {"OMP_NUM_THREADS": "1"}},
        "environment": {}, "checked_records": []}
    checks, launches = [], []
    def verify(*args):
        checks.append(args)
        return {} if outcome == "changed" and len(checks) > 1 else verified
    monkeypatch.setattr(driver, "verify_baseline", verify)
    def command(*args):
        argv = list(cell["argv"])
        argv[0] = "private-python"
        argv[argv.index("--output-directory") + 1] = str(tmp_path / driver.OUTPUT / "output")
        argv[argv.index("--json") + 1] = str(tmp_path / driver.OUTPUT / "metrics.json")
        argv.extend(["--checkpoint-source", str(baseline.parent)])
        return argv, []
    monkeypatch.setattr(driver, "control_command", command)
    monkeypatch.setattr(execution_module, "execution_environment", lambda *args: ({"LD_PRELOAD": "unsafe"}, {}))
    def lookup(*args):
        if outcome == "lookup":
            raise ValueError("Unexpected scientific import origin")
        return {"checked_records": []}
    monkeypatch.setattr(driver, "inspect_launcher", lookup)
    def execute(dataset, order, env, evidence, inputs, provenance):
        launches.append(dataset)
        assert Path.cwd() == launcher and order == ["private_phylogeny_control"]
        assert "LD_PRELOAD" not in env and env["PYTHONPATH"] == str(launcher)
        assert env["PYTHONNOUSERSITE"] == "1" and env["PYTHONDONTWRITEBYTECODE"] == "1"
        assert Path(env["NUMBA_CACHE_DIR"]).is_dir() and env["NUMBA_CACHE_LOCATOR_CLASSES"] == "UserProvidedCacheLocator"
        assert dataset["methods"][order[0]]["argv"][0] == "private-python"
        if outcome == "exception":
            raise RuntimeError("native interruption")
        new = tmp_path / driver.OUTPUT / "output/orthohmm_phylogeny"
        new.mkdir(parents=True)
        for name in driver.NATIVE_FILES:
            (new / name).write_text(name + ("different" if outcome == "parity" and name == driver.NATIVE_FILES[0] else ""))
        return {"failed_methods": order if outcome == "failed" else []}
    monkeypatch.setattr(execution_module, "execute", execute)
    cwd = Path.cwd()
    sha = "wrong" if outcome == "protocol" else record(protocol)["sha256"]
    if outcome == "success":
        driver.run(tmp_path, sha)
    else:
        with pytest.raises((ValueError, RuntimeError)):
            driver.run(tmp_path, sha)
    output = tmp_path / driver.OUTPUT
    if outcome == "protocol":
        assert not output.exists() and not checks and not launches
        return
    result = json.loads((output / "postflight.json").read_bytes())
    assert result["status"] == ("private_qfo_baseline_parity_complete_pending_admission" if outcome == "success" else "failed")
    assert all(result[key] is False for key in ("accuracy_evaluated", "publication_ready", "recovered_cpm_inference_authorized", "controlled_timing"))
    assert len(launches) == (0 if outcome == "lookup" else 1) and Path.cwd() == cwd
    assert cell == before and all((baseline / name).read_text() == name for name in driver.NATIVE_FILES)
    with pytest.raises(FileExistsError):
        driver.run(tmp_path, sha)


def test_no_scientific_import_before_gate():
    code = ("import sys; import benchmark_tools.run_qfo_private_phylogeny_control; "
            "assert not any(n == 'orthohmm' or n.startswith('orthohmm.') for n in sys.modules)")
    subprocess.run([sys.executable, "-S", "-B", "-c", code], check=True)


def test_pin_rejects_changed_file(tmp_path):
    path = tmp_path / "fixture.json"
    path.write_text("{}")
    value, pin = environment.read_pin(path, record(path)["sha256"])
    assert value == {} and pin == record(path)
    with pytest.raises(ValueError, match="identity"):
        environment.read_pin(path, "wrong")


@pytest.mark.parametrize("problem", [None, "process", "origin", "launcher", "shared_module", "shared_library"])
def test_launcher_probe_pins_actual_import_chain(tmp_path, monkeypatch, problem):
    launcher = tmp_path / "launcher"
    for path in (launcher / "orthohmm/phylogeny_pipeline.py", launcher / "benchmark_tools/replay_phylogeny.py", tmp_path / "private-python"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture")
    modules = {"orthohmm.phylogeny_pipeline": str(launcher / "orthohmm/phylogeny_pipeline.py"),
               "benchmark_tools.replay_phylogeny": str(launcher / "benchmark_tools/replay_phylogeny.py")}
    mapped = []
    if problem == "origin":
        modules["orthohmm.phylogeny_pipeline"] = str(tmp_path / "wrong-source.py")
    if problem == "launcher":
        modules["benchmark_tools.replay_phylogeny"] = str(tmp_path / "wrong-launcher.py")
    if problem == "shared_module":
        modules["unrelated"] = "/home/bizon/anaconda3/lib/python3.10/site-packages/unrelated.py"
    if problem == "shared_library":
        mapped = ["/home/bizon/.local/shared-library.so"]
    report = dict(modules=modules, mapped_files=mapped, executable=str(tmp_path / "private-python"))
    monkeypatch.setattr(environment.subprocess, "run", lambda command, **kwargs:
        subprocess.CompletedProcess(command, 1 if problem == "process" else 0, json.dumps(report), "diagnostic"))
    verified = {"launcher": str(launcher), "environment": {"tool_entrypoints": {
        "orthohmm_python": {"absolute_path": str(tmp_path / "private-python")}}}}
    output = tmp_path / "evidence"
    output.mkdir()
    if problem:
        with pytest.raises((ValueError, RuntimeError)):
            environment.inspect_launcher(verified, {}, output)
        assert (output / "lookup_process.json").is_file()
    else:
        result = environment.inspect_launcher(verified, {}, output)
        assert len(result["checked_records"]) == 3 and result["continuous_enforcement"] is False
        assert (output / "lookup.json").is_file()


def test_private_inventory_failure_propagates_before_environment_query(tmp_path, monkeypatch):
    from benchmark_tools import snapshot_runtime_trees as trees
    from benchmark_tools import run_simulation_methods as verifier

    legacy, ancestral, private, deployment = deployment_fixture()
    source = tmp_path / "record"
    source.write_text("fixture")
    pin = record(source)
    deployment.update(status="prospective_private_deployment_prepared", baseline=pin, previous_baseline=pin, candidate=pin, fixture=pin)
    candidate = {"status": "private_timing_environment_candidate_installed", "selected": {"numpy": "2.2.6"}}
    fixture = {"status": "both_native_prediction_fixtures_match", "candidate": pin}
    lookup = {"status": "native_lookup_repeated_identity_match", "baseline": pin}
    items = [deployment, private, ancestral, legacy, candidate, fixture, lookup, {"inventory": pin}, {}]
    def read(*args):
        return items.pop(0), pin
    monkeypatch.setattr(environment, "read_pin", read)
    def reject(*args):
        raise ValueError("Runtime inventory changed")
    monkeypatch.setattr(trees, "verify", reject)
    monkeypatch.setattr(verifier, "verify_environment", lambda *args: pytest.fail("must not query changed environment"))
    with pytest.raises(ValueError, match="Runtime inventory changed"):
        environment.verify_private_environment(tmp_path)


def test_launch_policy():
    path = Path(driver.__file__).parent / "results/qfo_private_phylogeny_control_20261001.sh"
    subprocess.run(["bash", "-n", str(path)], check=True)
    text = path.read_text()
    for token in ("--nodelist=bizon", "--cpus-per-task=32", "--mem=192G", "--time=24:00:00", "--no-requeue", "--protocol-sha256"):
        assert token in text
    assert "--dependency" not in text and "--array" not in text and "dgx" not in text
