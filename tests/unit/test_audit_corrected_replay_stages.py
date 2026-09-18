import copy
import json
import sys

import pytest

from benchmark_tools import audit_corrected_replay_stages as module
from benchmark_tools.checked_replay_payload_worker import FILES, STAGES
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


@pytest.mark.parametrize("problem", [None, "missing", "order", "execution", "manifest", "copy",
    "original_command", "command", "threads", "partition", "validation", "genes", "edges", "output"])
def test_stage_orchestration(tmp_path, monkeypatch, problem):
    # Real file hashes and stage handoff; native graph validation is mocked here
    # and covered with actual igraph/Leiden in test_validate_checked_replay_payload.
    root, executor, directory = tmp_path, tmp_path / "executor", tmp_path / "output"
    plan = tmp_path / "plan.json"
    write(plan, {})
    plan_record = record(plan)
    worker, observations = {"calls": []}, {}
    for index, label in enumerate(STAGES):
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
        manifest_path = stage / "payload_manifest.json"
        write(manifest_path, manifest)
        (stage / "partition.txt").write_text(f"gene_{index}\n")
        (stage / "worker.log").write_text("")
        partition = record(stage / "partition.txt")
        validation = {"status": "payload_checked", "accuracy_evaluated": False, "groups": 1,
                      "genes": 984137, "saved_graph": {"vertices": 984137, "edges": 10 + index},
                      "worker": {}, "provenance_checked": inputs}
        observations[str(payload)] = copy.deepcopy(validation)
        overrides = {k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
        row = {"index": index, "stage": label, "status": "checked", "exit_code": 0,
               "accuracy_evaluated": False, "manifest": record(manifest_path), "partition": partition,
               "validation": validation, "thread_environment": {"child_overrides": overrides,
                   "inherited": {**overrides, "OMP_NUM_THREADS": "32" if index >= 2 else "1"}},
               "command": [sys.executable, str(executor / "benchmark_tools/checked_replay_payload_worker.py"),
                   "--root", str(root), "--payload", str(payload), "--manifest", str(manifest_path),
                   "--manifest-sha256", record(manifest_path)["sha256"], "--corrected-plan", str(plan),
                   "--corrected-plan-sha256", plan_record["sha256"]]}
        write(stage / "execution.json", row)
        worker["calls"].append(row)
    replay = {"counts": {"rbnh_edges": 10, "multipass_edges": 11}, "stages": [
        {"label": label, "output": worker["calls"][index]["partition"]}
        for label, index in (("multipass", 1), ("multipass_refined", 1),
                             ("strict_profiles", 3), ("strict_profiles_refined", 3))]}
    def validate(payload, manifest, root, executor, corrected_plan, retained_partition):
        assert corrected_plan == plan_record
        assert retained_partition == record(payload.parent / "partition.txt")
        result = copy.deepcopy(observations[str(payload)])
        result["provenance_checked"].append(retained_partition)
        return result
    monkeypatch.setattr(module, "validate", validate)
    row = worker["calls"][0]
    stage = directory / "clustering/cluster_0_initial"
    if problem == "missing":
        worker["calls"].pop()
    elif problem == "order":
        worker["calls"].reverse()
    elif problem == "execution":
        write(stage / "execution.json", {})
    elif problem in ("manifest", "copy", "original_command"):
        manifest = json.loads((stage / "payload_manifest.json").read_text())
        if problem == "manifest":
            manifest["output_directory"] = "foreign"
        elif problem == "copy":
            manifest["original_payload"][0]["sha256"] = "changed"
        else:
            manifest["original_command"][-1] = "foreign"
        write(stage / "payload_manifest.json", manifest)
        row["manifest"] = record(stage / "payload_manifest.json")
    elif problem == "command":
        row["command"][-1] = "changed"
    elif problem == "threads":
        row["thread_environment"]["child_overrides"]["OMP_NUM_THREADS"] = "32"
    elif problem == "partition":
        (stage / "partition.txt").write_text("changed\n")
    elif problem == "validation":
        row["validation"]["groups"] = 2
    elif problem == "genes":
        row["validation"]["genes"] = 976504
        observations[str(stage / "payload")]["genes"] = 976504
    elif problem == "edges":
        replay["counts"]["rbnh_edges"] = 0
    elif problem == "output":
        replay["stages"][0]["output"] = row["partition"]
    if problem not in ("execution", "missing", "order"):
        write(stage / "execution.json", row)
    if problem:
        with pytest.raises(ValueError):
            module.audit(root, executor, directory, worker, replay, plan_record)
    else:
        result = module.audit(root, executor, directory, worker, replay, plan_record)
        assert result["status"] == "corrected_replay_stages_verified"
        assert len(result["clustering"]) == 4
        assert result["accuracy_evaluated"] is False
