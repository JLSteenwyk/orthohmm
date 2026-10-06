"""Retained admission launch evidence, not proof of live or scientific success."""

import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import run_assessment_gated_native_qfo_admission as gate
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


RESULTS = Path(gate.__file__).parent / "results"


@pytest.fixture
def receipts():
    names = dict(preflight="assessment_gated_native_qfo_admission_preflight_20261006.json",
        submission="assessment_gated_native_qfo_admission_submission_22452.json",
        release="assessment_gated_native_qfo_admission_released_22452.json",
        verified="assessment_gated_native_qfo_admission_verified_22452.json")
    return {k: json.loads((RESULTS / name).read_text()) for k, name in names.items()}


def test_actual_held_release_repoll_chain(receipts):
    sub, release, verified = (receipts[k] for k in ("submission", "release", "verified"))
    assert sub["held_inspection_passed"] is True and sub["released"] is False
    assert sub["preflight"] == record(RESULTS / "assessment_gated_native_qfo_admission_preflight_20261006.json")
    assert release["submission"] == verified["submission"] == record(
        RESULTS / "assessment_gated_native_qfo_admission_submission_22452.json")
    assert verified["initial_release"] == record(RESULTS / "assessment_gated_native_qfo_admission_released_22452.json")
    assert release["released"] is verified["released"] is True
    assert release["release"]["command"] == ["scontrol", "release", "22452"]
    assert release["release"]["returncode"] == 0


def test_launch_source_protocol_and_prepared_git_milestone_match(receipts):
    preflight, sub = receipts["preflight"], receipts["submission"]
    assert sub["source_commit"].startswith("1b0915a8")
    for name in ("worker", "validator", "batch", "protocol", "assessment_submission", "conversion_submission",
                 "reviewer_submission", "request", "plan", "scorer_environment", "retained_scorer_preflight"):
        assert sub[name] == preflight[name]
        check(sub[name])
    for name in ("worker", "batch", "protocol"):
        ref = sub[name]
        relative = str(Path(ref["path"]).relative_to(gate.ROOT))
        blob = subprocess.run(["git", "show", sub["source_commit"] + ":" + relative],
                              cwd=gate.ROOT, capture_output=True, check=True).stdout
        assert len(blob) == ref["bytes"] and hashlib.sha256(blob).hexdigest() == ref["sha256"]
    assert sub["validator"]["sha256"] == gate.ADMISSION_SHA


@pytest.mark.parametrize("role", ("preflight", "submission", "release", "verified"))
def test_no_scientific_admission_or_successor_claim_at_launch(receipts, role):
    receipt = receipts[role]
    for flag in ("validator_invoked", "accuracy_admitted", "automatic_retry", "next_identity_authorized", "publication_ready"):
        assert receipt[flag] is False
    if role == "preflight":
        assert receipt["future_output_read"] is False
    else:
        assert (receipt["job_id"], receipt["native_job_id"], receipt["reviewer_job_id"],
                receipt["conversion_job_id"], receipt["assessment_job_id"]) == (22452, 22444, 22445, 22450, 22451)


def test_all_five_job_identities_and_fixed_admission_limits(receipts):
    sub, release, verified = (receipts[k] for k in ("submission", "release", "verified"))
    held, initial, final = sub["held_fields"], release["post_release_fields"], verified["scheduler_observations"][0]["fields"]
    for fields in (held, initial, final):
        assert fields["JobId"] == "22452" and fields["JobState"] == "PENDING"
        assert fields["Dependency"] == "afterany:22451(unfulfilled)"
        assert fields["NumCPUs"] == fields["CPUs/Task"] == "2"
        assert fields["MinMemoryNode"] == "32G" and fields["TimeLimit"] == "06:00:00"
        assert fields["Requeue"] == fields["Restarts"] == "0"
        assert fields["Comment"] == sub["request"]["sha256"]
        assert fields["Command"] == sub["batch"]["path"]
    assert (held["Reason"], initial["Reason"], final["Reason"]) == ("JobHeldUser", "None", "Dependency")
    native, reviewer, conversion, assessment = [o["fields"] for o in verified["scheduler_observations"][1:]]
    assert native["JobId"] == "22444" and native["JobState"] == "RUNNING"
    for fields, job, predecessor in ((reviewer, 22445, 22444), (conversion, 22450, 22445), (assessment, 22451, 22450)):
        assert fields["JobId"] == str(job) and fields["JobState"] == "PENDING"
        assert fields["Dependency"] == f"afterany:{predecessor}(unfulfilled)"
    assert all(o["command"][:3] == ["scontrol", "show", "job"] for o in verified["scheduler_observations"])


def test_dated_environment_reuse_and_actual_runtime_are_not_conflated(receipts):
    preflight = receipts["preflight"]
    retained = json.loads(Path(preflight["retained_scorer_preflight"]["path"]).read_text())
    assert preflight["scorer_runtime_bytes_freshly_rechecked"] is False
    assert preflight["scorer_environment"] == retained["scorer_environment"]
    assert retained["scorer_runtime_records_checked"] == 704
    assert preflight["helper_sources_checked"] == 920
    runtime = preflight["runtime"]
    assert runtime["invocation"] != runtime["binary"]["path"]
    assert runtime["packages"] == dict(Bio="1.87", numpy="2.2.6", psutil="7.2.2")
    assert "unknown and potentially tool-dependent" in " ".join(preflight["limitations"])


def test_residual_original_report_is_not_sufficient_for_export(receipts):
    protocol = Path(receipts["submission"]["protocol"]["path"]).read_text()
    normalized = " ".join(protocol.split())
    assert "require this new job's successful terminal accounting and final accuracy-admitted wrapper gate" in normalized
    assert "residual original admitted report does not suffice" in normalized
    assert "not independent biological generalization, paired uncertainty or publication readiness" in normalized
