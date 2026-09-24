import copy
import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import audit_cpm_recovery_prefix as module
from benchmark_tools.checked_replay_payload_worker import FILES, STAGES
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


@pytest.fixture
def predecessor(tmp_path, monkeypatch):
    root, executor, directory = tmp_path, tmp_path / "executor", tmp_path / "output"
    plan = tmp_path / "plan.json"
    write(plan, {})
    plan_record = record(plan)
    worker = {"status": "failed", "accuracy_evaluated": False, "calls": []}
    observations = {}
    for index, label in enumerate(STAGES[:3]):
        stage = directory / "clustering" / f"cluster_{index}_{label}"
        payload = stage / "payload"
        payload.mkdir(parents=True)
        for name in FILES:
            (payload / name).write_text(name)
        inputs = [record(payload / name) for name in FILES]
        originals = [{**r, "path": str(tmp_path / "original" / name)} for name, r in zip(FILES, inputs)]
        manifest = {"stage": label, "index": index, "accuracy_evaluated": False, "inputs": inputs,
                    "original_payload": originals, "output_directory": str(directory / "replay"),
                    "original_command": [sys.executable, "-m", "orthohmm.leiden_worker", str(tmp_path / "original")]}
        write(stage / "payload_manifest.json", manifest)
        (stage / "partition.txt").write_text(f"gene_{index}\n")
        (stage / "worker.log").write_text("")
        validation = {"status": "payload_checked", "accuracy_evaluated": False, "groups": 1,
                      "genes": 984137, "saved_graph": {"vertices": 984137, "edges": 10 + index},
                      "worker": {}, "provenance_checked": inputs}
        observations[str(payload)] = copy.deepcopy(validation)
        manifest_record = record(stage / "payload_manifest.json")
        overrides = {k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
        row = {"index": index, "stage": label, "status": "checked", "exit_code": 0,
               "accuracy_evaluated": False, "manifest": manifest_record,
               "partition": record(stage / "partition.txt"), "validation": validation,
               "thread_environment": {"child_overrides": overrides,
                   "inherited": {**overrides, "OMP_NUM_THREADS": "32" if index == 2 else "1"}},
               "command": [sys.executable, str(executor / "benchmark_tools/checked_replay_payload_worker.py"),
                   "--root", str(root), "--payload", str(payload), "--manifest", manifest_record["path"],
                   "--manifest-sha256", manifest_record["sha256"], "--corrected-plan", str(plan),
                   "--corrected-plan-sha256", plan_record["sha256"], "--cpm-arm", "cpm_high"]}
        write(stage / "execution.json", row)
        worker["calls"].append(row)
    worker["calls"].append({"index": 3, "stage": STAGES[3], "status": "failed", "accuracy_evaluated": False})
    (directory / "replay").mkdir()
    (directory / "replay/orthogroups_multipass.txt").write_text("gene_1\n")

    def validate(payload, manifest, root, executor, **kwargs):
        partition = record(payload.parent / "partition.txt")
        assert kwargs == dict(corrected_plan=plan_record, retained_partition=partition, cpm_arm="cpm_high")
        result = copy.deepcopy(observations[str(payload)])
        result["provenance_checked"].append(partition)
        return result

    monkeypatch.setattr(module, "validate", validate)
    return root, executor, directory, worker, plan_record, observations


@pytest.mark.parametrize("problem", [None, "parent_success", "missing", "order", "failed_not_last",
    "boolean_index", "boolean_exit", "accuracy", "execution", "manifest", "copy", "payload",
    "original_command", "command", "threads", "partition", "validation", "genes", "output",
    "validator_failure", "mutation_during_validation"])
def test_predecessor_gates(predecessor, monkeypatch, problem):
    root, executor, directory, worker, plan_record, observations = predecessor
    row = worker["calls"][0]
    stage = directory / "clustering/cluster_0_initial"
    if problem == "parent_success":
        worker["status"] = "complete"
    elif problem == "missing":
        worker["calls"].pop()
    elif problem == "order":
        worker["calls"].reverse()
    elif problem == "failed_not_last":
        row["status"] = "failed"
    elif problem == "boolean_index":
        row["index"] = False
    elif problem == "boolean_exit":
        row["exit_code"] = False
    elif problem == "accuracy":
        worker["accuracy_evaluated"] = True
    elif problem == "execution":
        write(stage / "execution.json", {})
    elif problem in ("manifest", "copy", "original_command"):
        path = stage / "payload_manifest.json"
        manifest = json.loads(path.read_text())
        if problem == "manifest":
            manifest["output_directory"] = "foreign"
        elif problem == "copy":
            manifest["original_payload"][0]["sha256"] = "changed"
        else:
            manifest["original_command"][-1] = "foreign"
        write(path, manifest)
        row["manifest"] = record(path)
    elif problem == "payload":
        (stage / "payload/weights.npy").write_text("changed")
    elif problem == "command":
        row["command"][-1] = "cpm_low"
    elif problem == "threads":
        row["thread_environment"]["child_overrides"]["OMP_NUM_THREADS"] = "32"
    elif problem == "partition":
        (stage / "partition.txt").write_text("changed\n")
    elif problem == "validation":
        row["validation"]["groups"] = 2
    elif problem == "genes":
        row["validation"]["genes"] = 123
        observations[str(stage / "payload")]["genes"] = 123
    elif problem == "output":
        (directory / "replay/orthogroups_multipass.txt").write_text("wrong\n")
    elif problem == "validator_failure":
        def fail(*args, **kwargs):
            raise ValueError("native validation failed")
        monkeypatch.setattr(module, "validate", fail)
    elif problem == "mutation_during_validation":
        original = module.validate
        def mutate(*args, **kwargs):
            result = original(*args, **kwargs)
            (stage / "payload/weights.npy").write_text("changed after validation")
            return result
        monkeypatch.setattr(module, "validate", mutate)
    if problem != "execution":
        write(stage / "execution.json", row)
    if problem:
        with pytest.raises(ValueError):
            module.audit_prefix(root, executor, directory, worker, plan_record)
    else:
        result = module.audit_prefix(root, executor, directory, worker, plan_record)
        assert result["status"] == "cpm_recovery_predecessors_verified_not_admitted"
        assert len(result["clustering"]) == 3
        assert result["recovery_authorized"] is False
        assert result["accuracy_evaluated"] is False
        assert result["publication_ready"] is False
        assert worker["calls"][3]["status"] == "failed"


def test_run_refuses_existing_output(tmp_path):
    output = tmp_path / "report.json"
    output.write_text("preserve")
    with pytest.raises(FileExistsError):
        module.run(tmp_path, output)
    assert output.read_text() == "preserve"


def test_isolated_cli_help(tmp_path):
    done = subprocess.run([sys.executable, "-I", str(Path(module.__file__).resolve()), "--help"],
                          cwd=tmp_path, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert "--root" in done.stdout and "--output" in done.stdout
