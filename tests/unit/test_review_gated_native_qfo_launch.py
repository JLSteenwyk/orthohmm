"""Read back retained launch evidence, without claiming current scheduler state."""

import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import run_review_gated_native_qfo_pairs as gate
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


RESULTS = Path(gate.__file__).parent / "results"


def load(name):
    return json.loads((RESULTS / name).read_text())


@pytest.fixture
def receipts():
    return dict(
        preflight=load("review_gated_native_qfo_conversion_preflight_20261006.json"),
        submission=load("review_gated_native_qfo_conversion_submission_22450.json"),
        release=load("review_gated_native_qfo_conversion_released_22450.json"),
        verified=load("review_gated_native_qfo_conversion_verified_22450.json"))


def test_exact_launch_chain_and_no_repeated_release(receipts):
    sub, release, verified = (receipts[k] for k in ("submission", "release", "verified"))
    assert sub["preflight"] == record(RESULTS / "review_gated_native_qfo_conversion_preflight_20261006.json")
    assert release["submission"] == verified["submission"] == record(
        RESULTS / "review_gated_native_qfo_conversion_submission_22450.json")
    assert verified["initial_release"] == record(RESULTS / "review_gated_native_qfo_conversion_released_22450.json")
    assert sub["released"] is False and sub["held_inspection_passed"] is True
    assert release["released"] is verified["released"] is True
    assert release["release"]["command"] == ["scontrol", "release", "22450"]
    assert release["release"]["returncode"] == 0
    assert all(o["command"][:3] == ["scontrol", "show", "job"] for o in verified["scheduler_observations"])


def test_launch_pins_match_preflight_and_prepared_git_milestone(receipts):
    preflight, sub = receipts["preflight"], receipts["submission"]
    assert sub["source_commit"].startswith("b8192281")
    for name in ("worker", "batch", "protocol", "converter", "request", "plan", "reviewer_submission"):
        ref = sub[name]
        assert ref == preflight["source" if name == "worker" else name]
        check(ref)
    for name in ("worker", "batch", "protocol"):
        ref = sub[name]
        relative = str(Path(ref["path"]).relative_to(gate.ROOT))
        blob = subprocess.run(["git", "show", sub["source_commit"] + ":" + relative],
                              cwd=gate.ROOT, capture_output=True, check=True).stdout
        assert len(blob) == ref["bytes"] and hashlib.sha256(blob).hexdigest() == ref["sha256"]
    assert sub["converter"]["sha256"] == gate.CONVERTER_SHA


@pytest.mark.parametrize("role", ("preflight", "submission", "release", "verified"))
def test_launch_does_not_claim_scoring_successor_or_publication_completion(receipts, role):
    receipt = receipts[role]
    for name in ("accuracy_evaluated", "automatic_retry", "next_identity_authorized", "publication_ready"):
        assert receipt[name] is False
    if role == "preflight":
        assert receipt["future_review_read"] is receipt["conversion_started"] is False
    else:
        assert receipt["pair_conversion_completed"] is False
        assert (receipt["job_id"], receipt["native_job_id"], receipt["reviewer_job_id"]) == (22450, 22444, 22445)


def test_held_and_repolled_request_and_resource_envelopes_match(receipts):
    sub, verified = receipts["submission"], receipts["verified"]
    held, current = sub["held_fields"], verified["scheduler_observations"][0]["fields"]
    for fields in (held, current):
        assert fields["JobId"] == "22450" and fields["JobState"] == "PENDING"
        assert fields["Dependency"] == "afterany:22445(unfulfilled)"
        assert fields["NumCPUs"] == fields["CPUs/Task"] == "2"
        assert fields["MinMemoryNode"] == "32G" and fields["TimeLimit"] == "06:00:00"
        assert fields["Requeue"] == fields["Restarts"] == "0"
        assert fields["Comment"] == sub["request"]["sha256"]
        assert fields["Command"] == sub["batch"]["path"]
    assert held["Reason"] == "JobHeldUser" and current["Reason"] == "Dependency"
    native, reviewer = (o["fields"] for o in verified["scheduler_observations"][1:])
    assert native["JobId"] == "22444" and native["JobState"] == "RUNNING"
    assert reviewer["JobId"] == "22445" and reviewer["Dependency"] == "afterany:22444(unfulfilled)"


def test_initial_transient_status_is_retained_not_rewritten(receipts):
    release, verified = receipts["release"], receipts["verified"]
    fields = dict(token.split("=", 1) for token in release["post_release_controller"]["stdout"].split()
                  if "=" in token)
    assert fields["JobState"] == "PENDING" and fields["Reason"] == "None"
    assert fields["Dependency"] == "afterany:22445(unfulfilled)"
    assert "post_release_inspection_passed" not in release
    assert verified["post_release_inspection_passed"] is True
    assert "Re-poll original handle" in verified["initial_observation_assertion"]


def test_original_venv_not_resolved_binary_and_shared_capacity_are_explicit(receipts):
    for role in ("preflight", "submission", "release"):
        receipt = receipts[role]
        runtime = receipt["runtime"]
        assert runtime["invocation"] != runtime["binary"]["path"]
        assert runtime["prefix"] == str(Path(runtime["invocation"]).parent.parent)
        assert runtime["packages"] == dict(Bio="1.87", numpy="2.2.6", psutil="7.2.2")
        assert receipt["resource_context"]["available_memory_bytes"] >= 32 * 2**30
        assert receipt["resource_context"]["available_disk_bytes"] >= 128 * 2**30
    assert receipts["preflight"]["helper_sources_checked"] == receipts["submission"]["helper_sources_checked"] == 920
    assert "All system workloads" in receipts["preflight"]["host_cpu_snapshot"]["attribution"]
    assert "unknown and potentially tool-dependent" in " ".join(receipts["preflight"]["limitations"])
