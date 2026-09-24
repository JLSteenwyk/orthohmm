import json

import pytest

from benchmark_tools import prepare_cpm_checkpoint_recovery as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.mark.parametrize("problem", [None, "completed", "exit", "cpu", "node", "missing", "duplicate"])
def test_original_failure_scheduler(problem):
    fields = ["22081_1", "22081", "FAILED", "1:0", "32", "bizon"]
    for key, index, value in (("completed", 2, "COMPLETED"), ("exit", 3, "0:0"),
                             ("cpu", 4, "1"), ("node", 5, "other")):
        if key == problem:
            fields[index] = value
    rows = ["|".join(fields)]
    if problem == "missing":
        rows = []
    elif problem == "duplicate":
        rows *= 2
    text = "JobID|JobIDRaw|State|ExitCode|AllocCPUS|NodeList\n" + "\n".join(rows)
    if problem:
        with pytest.raises(ValueError):
            module.failed_job(text)
    else:
        assert module.failed_job(text)["JobIDRaw"] == "22081"


@pytest.mark.parametrize("problem", [None, "parent", "arm", "index", "job", "runtime", "source",
    "accuracy", "execution", "signal", "manifest", "stage", "bool_index"])
def test_failure_binding(tmp_path, problem):
    source = tmp_path / "benchmark_tools/run_qfo_cpm_variant.py"
    source.parent.mkdir()
    source.write_text("source")
    context = {"baseline_plan": {"runtime": {"verified": True}}}
    parent = dict(status="failed", arm="cpm_high", index=1, job_id="22081",
        executor_commit="a2486a9ef39afc2035c3fd20942696be48c5647d", context=context,
        runtime_before=context["baseline_plan"]["runtime"], source=record(source),
        accuracy_evaluated=False, publication_ready=False)
    manifest = {"path": "manifest", "sha256": "hash"}
    execution = dict(status="failed", index=3, stage="profile_expanded", error="SIGSEGV", manifest=manifest)
    worker = {"calls": [{}, {}, {}, dict(execution)]}
    if problem in ("parent", "arm", "index", "job", "runtime", "source", "accuracy", "bool_index"):
        key, value = {"parent": ("status", "complete"), "arm": ("arm", "cpm_low"), "index": ("index", 0),
            "job": ("job_id", "foreign"), "runtime": ("runtime_before", {}), "source": ("source", {}),
            "accuracy": ("accuracy_evaluated", True), "bool_index": ("index", True)}[problem]
        parent[key] = value
    elif problem == "execution":
        worker["calls"][3] = {}
    elif problem:
        key, value = {"signal": ("error", "other failure"), "manifest": ("manifest", {}),
                      "stage": ("stage", "multipass")}[problem]
        execution[key] = value
        worker["calls"][3] = dict(execution)
    args = (parent, execution, manifest, worker, context, {"JobIDRaw": "22081"}, tmp_path)
    if problem:
        with pytest.raises(ValueError):
            module.failure_identity(*args)
    else:
        module.failure_identity(*args)


def test_pending_refinement_creates_nothing(tmp_path, monkeypatch):
    def pending(root):
        raise ValueError("reconstruction pending")
    monkeypatch.setattr(module, "validate_refinement", pending)
    with pytest.raises(ValueError, match="pending"):
        module.prepare(tmp_path, tmp_path / "output")
    assert not (tmp_path / "output").exists()


def test_existing_output_is_never_reused(tmp_path, monkeypatch):
    def unexpected(root):
        raise AssertionError("Must reject before checking reconstruction")
    monkeypatch.setattr(module, "validate_refinement", unexpected)
    with pytest.raises(FileExistsError):
        module.prepare(tmp_path, tmp_path)


def test_predecessor_failure_is_retained(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "validate_refinement", lambda root: dict(
        status="original_cpm_refinement_independently_verified", comparison=dict(partition_equal=True, genes=984137),
        checked_records=[]))
    monkeypatch.setattr(module, "record", lambda path: dict(path=str(path), sha256=module.PROTOCOL_SHA, bytes=1))
    monkeypatch.setattr(module, "check", lambda item: None)
    monkeypatch.setattr(module, "HELPERS", {})
    monkeypatch.setattr(module, "ADAPTER_SHA", module.PROTOCOL_SHA)
    monkeypatch.setattr(module, "require_audit", lambda report: None)
    monkeypatch.setattr(module, "require_stop", lambda *args: None)
    def read(path, digest):
        if digest == module.AUDIT_SHA:
            return {"checked_records": []}
        if digest == module.BOUNDARY_SHA:
            return dict(status="frozen_worker_stopped_before_optimizer_unscored", job_id="22152", returncode=0,
                optimizer_called=False, observations=[dict(path="/fixture/constructor_adapter.json", sha256="adapter")],
                checked_records=[], worker_log=dict(path="/fixture/log", sha256="log"), result={"saved": {}})
        return {}
    monkeypatch.setattr(module, "read_frozen", read)
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobID|JobIDRaw|State|ExitCode|AllocCPUS|NodeList\n22081_1|22081|FAILED|1:0|32|bizon\n")
    def fail(root, output):
        raise ValueError("predecessor changed")
    monkeypatch.setattr(module, "audit_predecessors", fail)
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="predecessor changed"):
        module.prepare(tmp_path, output)
    status = json.loads((output / "status.json").read_text())
    assert status["status"] == "cpm_checkpoint_preflight_failed"
    assert status["preflight_passed"] is False
    assert status["optimizer_executed"] is False
