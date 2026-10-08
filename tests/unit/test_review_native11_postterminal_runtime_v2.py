import copy
from types import SimpleNamespace

import pytest

from benchmark_tools import native_factorial_allocated_execution as allocated
from benchmark_tools import review_native11_postterminal_runtime_v2 as current
from benchmark_tools import run_native_factorial_cost as historical


ROW = ("23985|COMPLETED|0:0|64|128G|bizon|gpu|1560|2026-10-07T15:27:53|"
       "2026-10-07T15:27:53|2026-10-08T08:56:00|orthohmm_allocated_factorial\n")


def test_verifier_is_original_allocation_aware_function_not_historical_route():
    assert current.verify_terminal is allocated.verify_terminal
    assert current.verify_terminal is not historical.verify_terminal
    fields = allocated.terminal_accounting(ROW, 23985)
    assert fields["AllocCPUS"] == "64"
    assert fields["JobName"] == "orthohmm_allocated_factorial"
    with pytest.raises(ValueError, match="Accounting identity"):
        historical.terminal_accounting(ROW, 23985)


def test_expired_controller_routes_through_original_allocated_accounting(monkeypatch):
    calls = []
    def run(command, **kwargs):
        calls.append(command)
        if command[0] == "scontrol":
            return SimpleNamespace(returncode=1, stdout="", stderr="Invalid job id specified")
        return SimpleNamespace(returncode=0, stdout=ROW, stderr="")
    monkeypatch.setattr(allocated.subprocess, "run", run)
    result = current.verify_terminal(23985)
    current.terminal_gate(result, {"sha256": "request"})
    assert result["source"] == "fresh_accounting_after_controller_expiry"
    assert [call[0] for call in calls] == ["scontrol", "sacct"]


@pytest.mark.parametrize("change", ["state", "exit", "comment"])
def test_success_and_live_request_comment_are_enforced(change):
    terminal = dict(source="live_controller",
        verified=dict(fields=dict(JobState="COMPLETED", ExitCode="0:0", Comment="request")))
    fields = terminal["verified"]["fields"]
    fields[{"state": "JobState", "exit": "ExitCode", "comment": "Comment"}[change]] = "wrong"
    with pytest.raises(ValueError):
        current.terminal_gate(terminal, {"sha256": "request"})


@pytest.mark.parametrize("change", ["job", "cpus", "memory", "name", "time"])
def test_original_allocation_parser_rejects_wrong_envelope(change):
    fields = ROW.strip().split("|")
    slot = {"job": 0, "cpus": 3, "memory": 4, "name": 11, "time": 7}[change]
    fields[slot] = "wrong"
    with pytest.raises(ValueError):
        allocated.terminal_accounting("|".join(fields), 23985)


def test_preflight_failure_is_retained_before_kernel_invocation(tmp_path, monkeypatch):
    destination = tmp_path / "component"
    monkeypatch.setattr(current, "DESTINATION", destination)
    refs = {}
    def record(path):
        key = str(path)
        if key == current.__file__:
            return dict(path=key, bytes=1, sha256="source")
        if key == current.kernel.__file__:
            return dict(path=key, bytes=1, sha256=current.KERNEL_SHA)
        digest = current.kernel.REQUEST_SHA if path == current.kernel.REQUEST else current.kernel.CLASSIFICATION_SHA
        return dict(path=key, bytes=1, sha256=digest)
    monkeypatch.setattr(current, "record", record)
    monkeypatch.setattr(current, "read", lambda ref: {})
    monkeypatch.setattr(current.kernel, "replay_runtime", lambda *a: pytest.fail("Must not replay"))
    def save(path, value):
        refs[path.name] = copy.deepcopy(value)
    monkeypatch.setattr(current, "save", save)
    with pytest.raises(ValueError, match="nonadmitting"):
        current.review(destination, "source")
    assert destination.is_dir()
    assert set(refs) == {"failure.json"}
    assert refs["failure.json"]["schema"] == current.SCHEMA
    for key in ("full_review_admitted", "next_identity_authorized", "automatic_retry", "publication_ready"):
        assert refs["failure.json"][key] is False


def test_existing_destination_refuses_before_source_reads(tmp_path, monkeypatch):
    destination = tmp_path / "component"
    destination.mkdir()
    monkeypatch.setattr(current, "DESTINATION", destination)
    monkeypatch.setattr(current, "record", lambda *a: pytest.fail("Should not read source"))
    with pytest.raises(ValueError, match="fresh fixed"):
        current.review(destination, "source")


def test_wrong_source_digest_refuses_before_output_creation(tmp_path, monkeypatch):
    destination = tmp_path / "component"
    monkeypatch.setattr(current, "DESTINATION", destination)
    monkeypatch.setattr(current, "record", lambda *a: {"sha256": "wrong"})
    with pytest.raises(ValueError, match="sources changed"):
        current.review(destination, "source")
    assert not destination.exists()
