import json
from pathlib import Path

import pytest

from benchmark_tools import finish_native_orthobench_attempt as followup


@pytest.fixture
def sample(tmp_path, monkeypatch):
    monkeypatch.setattr(followup, "ROOT", tmp_path)
    request = {"plan": {"path": "plan"}, "job_id": 22431, "index": 4}
    request_ref = {"path": "request"}
    events = []
    terminal = {"status": "native_success"}
    terminal_ref = {"path": "terminal"}
    def read(ref):
        return request if ref == request_ref else terminal if ref == terminal_ref else {}
    monkeypatch.setattr(followup, "read", read)
    monkeypatch.setattr(followup, "validate_request", lambda *args: events.append("validate_request"))
    monkeypatch.setattr(followup, "validate_plan", lambda *args: [{}] * 4 + [{"dataset": "orthobench"}])
    def review(ref, path):
        events.append("terminal_review")
        path.mkdir()
        return terminal_ref
    def score(ref, reviewed, path):
        assert ref == request_ref and reviewed == terminal_ref
        events.append("separate_score")
        path.mkdir()
        return {"path": "score"}
    monkeypatch.setattr(followup, "review", review)
    monkeypatch.setattr(followup, "score_attempt", score)
    monkeypatch.setattr(followup, "check", lambda *args: events.append("source_check"))
    return tmp_path, request_ref, terminal, events


def test_review_precedes_separate_score_without_inference_or_successor(sample):
    root, request, _, events = sample
    result = followup.finish(request, root / "review", root / "score")
    assert events == ["validate_request", "terminal_review", "separate_score", "source_check"]
    assert result["status"] == "terminal_reviewed_and_scored"
    assert result["parent_native_job_id"] == 22431
    assert not result["inference_reexecuted"] and not result["next_identity_released"]
    assert not result["automatic_retry"] and not result["publication_ready"]
    assert json.loads((root / "review/followup.json").read_text()) == result


def test_failed_native_review_is_retained_without_scoring(sample):
    root, request, terminal, events = sample
    terminal["status"] = "native_failure_retained"
    result = followup.finish(request, root / "review", root / "score")
    assert "separate_score" not in events and not (root / "score").exists()
    assert result["score"] is None and result["status"] == "native_failure_retained_unscored"


def test_unknown_outcome_is_not_relabelled_as_a_native_failure(sample):
    root, request, terminal, events = sample
    terminal["status"] = "unknown"
    with pytest.raises(ValueError, match="Unknown"):
        followup.finish(request, root / "review", root / "score")
    assert "separate_score" not in events


def test_live_or_failed_review_never_calls_scorer(sample, monkeypatch):
    root, request, _, events = sample
    def refuse(*args):
        events.append("review_refused")
        raise ValueError("not terminal")
    monkeypatch.setattr(followup, "review", refuse)
    with pytest.raises(ValueError, match="not terminal"):
        followup.finish(request, root / "review", root / "score")
    assert "separate_score" not in events and not (root / "score").exists()


@pytest.mark.parametrize("kind", ["same", "existing", "outside", "parent_alias", "symlink"])
def test_unsafe_or_used_destinations_refused_before_review(sample, kind):
    root, request, _, events = sample
    review, score = root / "review", root / "score"
    if kind == "same":
        score = review
    elif kind == "existing":
        score.mkdir()
    elif kind == "outside":
        score = root.parent / "outside"
    elif kind == "parent_alias":
        score = root / "child/../score"
    else:
        score.symlink_to(root / "absent")
    with pytest.raises(ValueError):
        followup.finish(request, review, score)
    assert "terminal_review" not in events


def test_qfo_is_not_silently_scored_as_orthobench(sample, monkeypatch):
    root, request, _, events = sample
    monkeypatch.setattr(followup, "validate_plan", lambda *args: [{}] * 4 + [{"dataset": "qfo_corrected"}])
    with pytest.raises(ValueError, match="QfO"):
        followup.finish(request, root / "review", root / "score")
    assert "terminal_review" not in events


def test_batch_envelope_matches_noncomparative_postprocessing_scope():
    script = Path(followup.__file__).with_suffix(".sh").read_text()
    assert "--cpus-per-task=2" in script and "--mem=8G" in script
    assert "--time=02:00:00" in script and "--no-requeue" in script
    assert "native_factorial_review_py310_20261004" in script
    assert "sbatch" not in script and "scontrol" not in script


@pytest.mark.parametrize("kind", ["wrong_python", "wrong_source", "wrong_request"])
def test_cli_refuses_incompatible_runtime_or_changed_pins(tmp_path, monkeypatch, kind):
    request = tmp_path / "request.json"
    request.write_text("{}")
    source_sha = followup.record(followup.__file__)["sha256"]
    request_sha = followup.record(request)["sha256"]
    argv = ["finish", "--request", str(request), "--request-sha256", request_sha,
            "--source-sha256", source_sha, "--review-directory", str(tmp_path / "review"),
            "--score-directory", str(tmp_path / "score")]
    monkeypatch.setattr(followup.sys, "version_info", (3, 12) if kind == "wrong_python" else (3, 10))
    if kind == "wrong_source":
        argv[6] = "0" * 64
    elif kind == "wrong_request":
        argv[4] = "0" * 64
    monkeypatch.setattr(followup.sys, "argv", argv)
    with pytest.raises(ValueError):
        followup.main()
    assert not (tmp_path / "review").exists() and not (tmp_path / "score").exists()
