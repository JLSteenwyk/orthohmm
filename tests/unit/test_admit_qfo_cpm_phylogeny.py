from copy import deepcopy
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import admit_qfo_cpm_phylogeny as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_admit_qfo_parameter_phylogeny import execution_fixture


def accounting(state="COMPLETED", exit_code="0:0", node="bizon", cpus="32"):
    return ("JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n"
            f"22088_0|23000|{state}|{exit_code}|01:00:00|{node}|{cpus}\n")


@pytest.mark.parametrize("index", [-1, 2, True, "0", 0.0])
def test_invalid_index(index):
    with pytest.raises(ValueError):
        module.completed_task(accounting(), index)


@pytest.mark.parametrize("change", [{"state": "RUNNING"}, {"state": "FAILED"},
    {"state": "PENDING"}, {"exit_code": "1:0"}, {"node": "spark-7ff0"}, {"cpus": "2"}])
def test_wrong_accounting_rejected(change):
    with pytest.raises(ValueError, match="completed"):
        module.completed_task(accounting(**change), 0)


def test_exact_unique_completed_task():
    assert module.completed_task(accounting(), 0)["JobIDRaw"] == "23000"
    for value, index in ((accounting() + accounting().splitlines()[1] + "\n", 0), (accounting(), 1)):
        with pytest.raises(ValueError):
            module.completed_task(value, index)


def test_live_job_precedes_artifact_reads(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: accounting(state="RUNNING"))
    monkeypatch.setattr(module, "verify_sources", lambda *a: pytest.fail("Read unfinished output"))
    with pytest.raises(ValueError, match="completed"):
        module.admit(tmp_path, 0)


@pytest.mark.parametrize("problem", [None, "extra", "missing", "command", "fresh", "status", "provenance"])
def test_postflight_requires_exact_fresh_admission(problem):
    status, pre, post, expected, manifest, cell = execution_fixture()
    command, fresh = ["python", "admitter.py"], {"path": "fresh.json", "sha256": "abc", "bytes": 1}
    post.update(admission_command=deepcopy(command), fresh_candidate_admission=deepcopy(fresh))
    if problem == "extra":
        post["unexpected"] = True
    elif problem == "missing":
        del post["fresh_candidate_admission"]
    elif problem == "command":
        post["admission_command"] = ["other.py"]
    elif problem == "fresh":
        post["fresh_candidate_admission"]["sha256"] = "wrong"
    elif problem == "status":
        status["failed_methods"] = [cell["label"]]
    elif problem == "provenance":
        status["provenance"]["job_id"] = "other"
    if problem:
        with pytest.raises(ValueError):
            module.check_execution(status, pre, post, expected, manifest, cell, command, fresh)
    else:
        module.check_execution(status, pre, post, expected, manifest, cell, command, fresh)


@pytest.mark.parametrize("problem", [None, "seed_partition", "candidate_partition", "membership_constraints"])
def test_adaptation_uses_complete_cpm_arm_without_mutation(problem):
    candidate = {"seed_partition": {"path": "CPM_seed"}, "candidate_partition": {"path": "CPM_candidate"},
        "membership_constraints": {"path": "CPM_constraints"}, "expansion": {"merges": 5}, "candidate_expansion": True}
    manifest = {"input_fastas": [], "candidate_arms": {"p1_c1": {"seed_partition": {"path": "baseline_seed"}},
                                                      "p0_c0": {"candidate_partition": {"path": "baseline"}}}}
    verified = {"manifest": manifest, "admission": {"candidate_arm": candidate}, "arm": {
        "seed_partition": deepcopy(candidate["seed_partition"]), "partition": deepcopy(candidate["candidate_partition"]),
        "constraints": deepcopy(candidate["membership_constraints"])}}
    if problem:
        candidate[problem] = {"path": "wrong"}
        with pytest.raises(ValueError):
            module.adapted_manifest(verified, {"path": "launcher"})
    else:
        before = deepcopy(verified)
        adapted = module.adapted_manifest(verified, {"path": "launcher"})
        assert adapted["candidate_arms"]["p1_c1"] == candidate
        assert adapted["candidate_arms"]["p0_c0"] == manifest["candidate_arms"]["p0_c0"]
        adapted["candidate_arms"]["p1_c1"]["seed_partition"]["path"] = "changed"
        assert verified == before


@pytest.mark.parametrize("problem", [None, "revision", "source", "helpers", "fresh", "pairs", "changed"])
def test_admission_orchestration(tmp_path, monkeypatch, problem):
    def write(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))
        return record(path)
    executor = tmp_path / "benchmarks/work/publication_qfo_cpm_phylogeny_v3"
    producer = write(executor / "benchmark_tools/run_qfo_cpm_phylogeny.py", "source")
    helper = write(executor / "benchmark_tools/helper.py", "helper")
    monkeypatch.setattr(module, "EXECUTOR_SHA", "wrong" if problem == "source" else producer["sha256"])
    output = tmp_path / "benchmarks/results/qfo_cpm_phylogeny_v1/cpm_low"
    candidate = {k: write(tmp_path / f"{k}.txt", k) for k in (
        "seed_partition", "candidate_partition", "membership_constraints")}
    verified = {"manifest": {"input_fastas": [], "candidate_arms": {}}, "original": {},
        "arm": {"seed_partition": candidate["seed_partition"], "partition": candidate["candidate_partition"],
                "constraints": candidate["membership_constraints"]},
        "admission": {"candidate_arm": candidate}, "admission_executor": str(tmp_path / "admitter"),
        "environment": {}, "context": {"arm": "cpm_low"}, "checked_records": list(candidate.values()),
        "launcher": str(tmp_path / "launcher"), "prepared": str(tmp_path / "prepared")}
    verified["admission_record"] = write(tmp_path / "admission.json", verified["admission"])
    checks = []
    def verify(*args):
        checks.append(args)
        return {} if problem == "changed" and len(checks) > 1 else verified
    monkeypatch.setattr(module, "verify_sources", verify)
    monkeypatch.setattr(module.subprocess, "check_output", lambda command, **k: accounting() if command[0] == "sacct"
        else "wrong" if problem == "revision" else module.EXECUTOR_COMMIT)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    cell = {"label": "candidate_cpm_low"}
    argv = [sys.executable, "launcher.py"]
    equivalence = [{"prepared": write(tmp_path / "prepared.py", "launcher"),
                    "executed": write(tmp_path / "executed.py", "launcher")}]
    monkeypatch.setattr(module, "variant_cell", lambda *a: cell)
    monkeypatch.setattr(module, "native_command", lambda *a: (argv, equivalence))
    monkeypatch.setattr(module, "execution_environment", lambda *a: ({}, {}))
    pre = {"source": producer, "helpers": [helper, producer], "verified": verified, "cell": cell,
        "executed_argv": argv, "launcher_source_equivalence": equivalence, "resolved_tools": {},
        "cwd": verified["launcher"], "job_id": "23000", "array_task_id": "0",
        "scope": "CPM-specific seed/candidate families; independently inferred phylogeny; no accuracy scoring"}
    if problem == "helpers":
        pre["helpers"].pop()
    fresh = write(output / "fresh_candidate_admission.json", {} if problem == "fresh" else verified["admission"])
    command = [sys.executable, "-B", str(tmp_path / "admitter/benchmark_tools/admit_qfo_cpm_candidates.py"),
        "--root", str(tmp_path), "--index", "0", "--output", fresh["path"]]
    post = {"status": "complete_pending_native_validation", "cell": cell, "accuracy_evaluated": False,
        "native_outputs_validated": False, "admission_command": command, "fresh_candidate_admission": fresh}
    status = {"provenance": pre, "verified_inputs": {"status": "ready", "inputs": []},
        "dataset": cell["label"], "methods": {cell["label"]: {}}, "status": "finished_pending_native_validation",
        "failed_methods": [], "accuracy_evaluated": False, "native_outputs_validated": False}
    for name, value in (("preflight.json", pre), ("postflight.json", post), ("execution/status.json", status),
        ("admission.log", "completed"), ("output/orthohmm_phylogeny/provenance_manifest.json", {"results": {"ortholog_pairs": 1}}),
        ("output/orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv", "pairs")):
        write(output / name, value)
    calls = []
    monkeypatch.setattr(module, "verify_process", lambda *a: calls.append("process") or [{}])
    def native(prepared, *args, **kwargs):
        calls.append("native")
        assert prepared["candidate_arms"]["p1_c1"] == candidate
        assert kwargs["expected_revision"] == module.LAUNCHER_COMMIT
        return {"status": "native_group_output_verified"}
    monkeypatch.setattr(module, "validate_native_cell", native)
    monkeypatch.setattr(module, "gene_ownership", lambda manifest, native, path: ({}, {}))
    monkeypatch.setattr(module, "check_pairs", lambda *a: 0 if problem == "pairs" else 1)
    if problem:
        with pytest.raises(ValueError):
            module.admit(tmp_path, 0)
    else:
        result = module.admit(tmp_path, 0)
        assert result["status"] == "cpm_native_pairs_verified_unscored"
        assert result["native_pair_count"] == 1 and result["context"] == verified["context"]
        assert result["accuracy_evaluated"] is result["scoring_admitted"] is result["publication_ready"] is False
        assert calls == ["process", "native", "process"] and len(checks) == 2
