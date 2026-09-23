import copy
import json
import sys

import pytest

from benchmark_tools import admit_qfo_cpm_candidates as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.mark.parametrize("problem", [None, "missing", "duplicate", "running", "exit", "node", "cpus", "index"])
def test_scheduler(problem):
    row = ["22084_0", "22085", "COMPLETED", "0:0", "00:10:00", "bizon", "2"]
    if problem == "running":
        row[2] = "RUNNING"
    elif problem == "exit":
        row[3] = "1:0"
    elif problem == "node":
        row[5] = "other"
    elif problem == "cpus":
        row[6] = "32"
    text = "JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n"
    text += ("|".join(row) + "\n") * (0 if problem == "missing" else 2 if problem == "duplicate" else 1)
    if problem:
        with pytest.raises(ValueError):
            module.completed_task(text, True if problem == "index" else 0)
    else:
        assert module.completed_task(text, 0)["JobIDRaw"] == "22085"


@pytest.mark.parametrize("index", [0, 1])
@pytest.mark.parametrize("problem", [None, "status", "arm", "index", "job", "executor", "source", "inputs", "context",
    "accuracy", "label", "delta", "parameters", "calls", "report", "expansion", "nan", "negative"])
def test_report(index, problem):
    params, source, inputs, context = {"min_norm": .03, "min_margin": 1.5}, {"source": True}, [{"input": True}], {"arm": index}
    expansion = {"parameters": copy.deepcopy(params)}
    report = {"status": "cpm_candidates_prepared_pending_admission", "arm": module.ARMS[index], "index": index,
        "array_task_id": str(index), "job_id": "fixture", "executor_commit": module.EXECUTOR, "source": source,
        "inputs": inputs, "context": context, "accuracy_evaluated": False, "publication_ready": False,
        "candidate_parameter_control": {"label": "control", "delta": {}, "applied_parameters": copy.deepcopy(params),
            "engine_calls": 1, "engine_fixed_profile_report": copy.deepcopy(expansion)},
        "candidate_arm": {"expansion": copy.deepcopy(expansion)}, "incremental_seconds": 1.0}
    control = report["candidate_parameter_control"]
    mutations = {"status": (report, "status", "failed"), "arm": (report, "arm", "control"),
        "index": (report, "index", bool(index)), "job": (report, "job_id", "other"),
        "executor": (report, "executor_commit", "other"), "source": (report, "source", {}),
        "inputs": (report, "inputs", []), "context": (report, "context", {}),
        "accuracy": (report, "accuracy_evaluated", True), "label": (control, "label", "norm_low"),
        "delta": (control, "delta", {"min_norm": .024}), "parameters": (control, "applied_parameters", {}),
        "calls": (control, "engine_calls", True), "report": (control, "engine_fixed_profile_report", {"parameters": {}}),
        "expansion": (report["candidate_arm"], "expansion", {}), "nan": (report, "incremental_seconds", float("nan")),
        "negative": (report, "incremental_seconds", -1)}
    args = (report, index, {"JobIDRaw": "fixture"}, context, params, source, inputs)
    if problem:
        target, key, value = mutations[problem]
        target[key] = value
        with pytest.raises(ValueError):
            module.validate_report(*args)
    else:
        module.validate_report(*args)


@pytest.mark.parametrize("state", ["PENDING", "RUNNING", "FAILED"])
def test_active_or_failed_job_precedes_file_access(tmp_path, monkeypatch, state):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        f"JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n22084_0|22085|{state}|0:0|00:01:00|bizon|2\n")
    output = tmp_path / "admission.json"
    with pytest.raises(ValueError, match="completed"):
        module.admit(tmp_path, 0, output)
    assert not output.exists()


@pytest.mark.parametrize("problem", [None, "fresh", "runtime", "numeric", "auditor", "content", "merge", "source"])
def test_orchestration(tmp_path, monkeypatch, problem):
    def write(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))
        return record(path)
    root = tmp_path
    executor = root / "benchmarks/work/publication_qfo_cpm_candidates_v4"
    source = write(executor / "benchmark_tools/prepare_qfo_cpm_candidates.py", {})
    auditor = write(executor / "benchmark_tools/audit_accuracy_checkpoint.py", {})
    monkeypatch.setattr(module, "SOURCE_SHA", "wrong" if problem == "source" else source["sha256"])
    directory = root / "benchmarks/results/qfo_cpm_candidates_v1/cpm_low"
    names = root / "names.txt"
    names.write_text("".join(f"g{i}\n" for i in range(984137)))
    checkpoint = write(root / "checkpoint/manifest.json", {})
    plan = {"runtime": {}, "checkpoint_manifest": checkpoint}
    plan_record = write(root / "plan.json", plan)
    numeric = {"status": "numeric_checkpoint_verified", "manifest": checkpoint,
               "summary": {"genes": 984137, "species": 78}, "accuracy_evaluated": False, "auditor": auditor}
    baseline = {"runtime_before": {}, "runtime_after": {}, "numeric_checkpoint": numeric,
        "candidate_arms": {"p1_c1": {"expansion": {"parameters": {"min_norm": .03, "min_margin": 1.5}}}}}
    baseline_record = write(root / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json", baseline)
    monkeypatch.setattr(module, "BASELINE_SHA", baseline_record["sha256"])
    write(root / "benchmark_tools/results/publication_native_runtime_20260916.json", {})
    context = {"cwd": str(root / "launcher"), "checked_records": []}
    seed = write(root / "seed.txt", {})
    candidate = write(directory / "candidate/partition.txt", {})
    trace = write(directory / "candidate/merges.json", [])
    arm = {"seed_partition": seed, "candidate_partition": candidate, "membership_constraints": trace,
           "output_files": [candidate, trace]}
    replay = {"executor": str(root / "replay_validator"), "report": {"fixture": True},
              "report_record": write(root / "replay_admission.json", {}), "seed_partition": seed, "checked_records": [seed]}
    fresh = write(directory / "fresh_replay_admission.json", {} if problem == "fresh" else replay["report"])
    content = {"checked_records": [seed, candidate, trace], "status": "fixture"}
    prior_numeric = copy.deepcopy(numeric)
    if problem == "numeric":
        prior_numeric["summary"] = {}
    elif problem == "auditor":
        prior_numeric["auditor"] = {}
    report = {"replay_admission": replay, "fresh_replay_admission": fresh,
        "admission_command": [sys.executable, "-B", str(root / "replay_validator/benchmark_tools/admit_qfo_cpm_variant.py"),
                              "--root", str(root), "--index", "0", "--output", fresh["path"]],
        "runtime_before": {}, "runtime_after": {}, "numeric_checkpoint": prior_numeric,
        "candidate_arm": arm, "content_audit": {} if problem == "content" else content}
    write(directory / "manifest.json", report)
    write(directory / "admission.log", {})
    write(root / "benchmarks/work/qfo_cpm_candidates_22084_0.time.txt", {})
    monkeypatch.setattr(module, "corrected_evidence", lambda *a: (plan, plan_record, {}, record(names)))
    monkeypatch.setattr(module, "evidence", lambda *a: context)
    monkeypatch.setattr(module, "replay_evidence", lambda *a: replay)
    monkeypatch.setattr(module, "validate_report", lambda *a: None)
    monkeypatch.setattr(module, "verify", lambda *a: {"changed": True} if problem == "runtime" else {})
    monkeypatch.setattr(module, "audit_numeric", lambda *a: numeric)
    def audit(*args):
        assert args[0] == arm and args[1] == seed and len(args[3]) == 984137 and args[4] is True
        return content
    monkeypatch.setattr(module, "audit_arm", audit)
    monkeypatch.setattr(module, "partition", lambda *a: ([set("ab")], {}))
    reconstructions = []
    def reconstruct(*args):
        reconstructions.append(args)
        if problem == "merge":
            raise ValueError("Invalid merge reconstruction")
    monkeypatch.setattr(module, "validate_merge_reconstruction", reconstruct)
    def check_output(command, **kwargs):
        if command[0] == "sacct":
            return "JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n22084_0|22085|COMPLETED|0:0|00:01:00|bizon|2\n"
        return module.EXECUTOR
    monkeypatch.setattr(module.subprocess, "check_output", check_output)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    output = root / "admitted.json"
    if problem:
        with pytest.raises(ValueError):
            module.admit(root, 0, output)
        assert not output.exists()
    else:
        result = module.admit(root, 0, output)
        assert json.loads(output.read_text()) == result
        assert result["status"] == "cpm_candidates_admitted_unscored"
        assert result["candidate_arm"] == arm and len(reconstructions) == 1
        assert result["accuracy_evaluated"] is False and result["publication_ready"] is False
        with pytest.raises(FileExistsError):
            module.admit(root, 0, output)
