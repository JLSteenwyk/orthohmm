import json

import pytest

from benchmark_tools.audit_threadripper_native_outcome import audit


def setup(tmp_path, outcome="exited_zero"):
    command = ["/native/python", "-m", "orthohmm"]
    done = dict(exit_code=0 if outcome == "exited_zero" else 1,
                timed_out=outcome == "timed_out")
    (tmp_path / "done.json").write_text(json.dumps(done))
    (tmp_path / "command.json").write_text(json.dumps(dict(command=command)))
    (tmp_path / "native.log").write_text("retained output")
    run = dict(measurement_directory=str(tmp_path), native_argv=command, cwd="/frozen")
    replay = dict(measured=dict(native=done), native_outcome=outcome)
    return run, replay


def test_success_replays_expected_command_and_checks_outputs(tmp_path):
    run, replay = setup(tmp_path)
    calls = []
    def reproduce(path, job, command):
        assert path == tmp_path and job == 42 and command == run["native_argv"]
        calls.append("replay")
        return replay
    def outputs(actual, measured, baseline):
        assert actual == run and measured["cwd"] == "/frozen" and baseline == {"frozen": True}
        calls.append("outputs")
        return dict(status="threadripper_native_outputs_checked")
    result = audit(run, {"frozen": True}, 42, replay_fn=reproduce, validate_fn=outputs)
    assert calls == ["replay", "outputs"]
    assert result["status"] == "native_success_outputs_verified"
    assert result["next_submission_authorized"] is False


@pytest.mark.parametrize("outcome", ["exited_nonzero", "timed_out"])
def test_failure_is_retained_without_output_validation_or_retry(tmp_path, outcome):
    run, replay = setup(tmp_path, outcome)
    def no_outputs(*args):
        pytest.fail("Must not require successful outputs for native failure")
    result = audit(run, {}, 42, replay_fn=lambda *a: replay, validate_fn=no_outputs)
    assert result["status"] == "native_failure_requires_review"
    assert result["native_outcome"] == outcome
    assert result["automatic_retry"] is False


def test_changed_command_prevents_replay(tmp_path):
    run, _ = setup(tmp_path)
    run["native_argv"] = ["/other"]
    with pytest.raises(ValueError, match="command"):
        audit(run, {}, 42, replay_fn=lambda *a: pytest.fail("Unexpected replay"))


def test_replay_mismatch_rejected(tmp_path):
    run, replay = setup(tmp_path)
    replay["measured"]["native"]["exit_code"] = 1
    with pytest.raises(ValueError, match="outcome differs"):
        audit(run, {}, 42, replay_fn=lambda *a: replay)


def test_evidence_mutation_during_validation_rejected(tmp_path):
    run, replay = setup(tmp_path)
    def outputs(*args):
        (tmp_path / "native.log").write_text("changed")
        return dict(status="threadripper_native_outputs_checked")
    with pytest.raises(ValueError, match="identity changed"):
        audit(run, {}, 42, replay_fn=lambda *a: replay, validate_fn=outputs)


def test_output_validation_failure_propagates(tmp_path):
    run, replay = setup(tmp_path)
    def outputs(*args):
        raise ValueError("Missing native groups")
    with pytest.raises(ValueError, match="Missing native"):
        audit(run, {}, 42, replay_fn=lambda *a: replay, validate_fn=outputs)
