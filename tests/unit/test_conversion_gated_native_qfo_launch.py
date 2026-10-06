"""Verify retained assessment launch evidence, not current or future execution."""

import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import run_conversion_gated_native_qfo_assessment as gate
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


RESULTS = Path(gate.__file__).parent / "results"


@pytest.fixture
def receipts():
    names = dict(preflight="conversion_gated_native_qfo_assessment_preflight_20261006.json",
        submission="conversion_gated_native_qfo_assessment_submission_22451.json",
        release="conversion_gated_native_qfo_assessment_released_22451.json",
        verified="conversion_gated_native_qfo_assessment_verified_22451.json")
    return {key: json.loads((RESULTS / name).read_text()) for key, name in names.items()}


def test_actual_held_release_and_repoll_chain(receipts):
    sub, release, verified = (receipts[k] for k in ("submission", "release", "verified"))
    assert sub["held_inspection_passed"] is True and sub["released"] is False
    assert sub["preflight"] == record(RESULTS / "conversion_gated_native_qfo_assessment_preflight_20261006.json")
    assert release["submission"] == verified["submission"] == record(
        RESULTS / "conversion_gated_native_qfo_assessment_submission_22451.json")
    assert verified["initial_release"] == record(RESULTS / "conversion_gated_native_qfo_assessment_released_22451.json")
    assert release["released"] is verified["released"] is True
    assert release["release"]["command"] == ["scontrol", "release", "22451"]
    assert release["release"]["returncode"] == 0


def test_current_and_committed_source_pins_match_launch(receipts):
    preflight, sub = receipts["preflight"], receipts["submission"]
    assert sub["source_commit"].startswith("4d532428")
    for name in ("worker", "driver", "batch", "protocol", "conversion_submission", "reviewer_submission",
                 "request", "plan", "scorer_environment"):
        assert sub[name] == preflight[name]
        check(sub[name])
    for name in ("worker", "batch", "protocol"):
        ref = sub[name]
        relative = str(Path(ref["path"]).relative_to(gate.ROOT))
        blob = subprocess.run(["git", "show", sub["source_commit"] + ":" + relative],
                              cwd=gate.ROOT, capture_output=True, check=True).stdout
        assert len(blob) == ref["bytes"] and hashlib.sha256(blob).hexdigest() == ref["sha256"]
    assert sub["driver"]["sha256"] == gate.ASSESSMENT_SHA


@pytest.mark.parametrize("role", ("preflight", "submission", "release", "verified"))
def test_no_completed_assessment_or_admission_claim(receipts, role):
    receipt = receipts[role]
    for name in ("accuracy_admitted", "automatic_retry", "next_identity_authorized", "publication_ready"):
        assert receipt[name] is False
    if role == "preflight":
        assert receipt["future_conversion_read"] is receipt["assessment_driver_invoked"] is False
        assert receipt["command_executed"] is False
    else:
        assert receipt["assessment_completed"] is False
        assert (receipt["job_id"], receipt["native_job_id"], receipt["reviewer_job_id"],
                receipt["conversion_job_id"]) == (22451, 22444, 22445, 22450)


def test_actual_all_four_job_identities_and_fixed_assessment_limits(receipts):
    sub, release, verified = (receipts[k] for k in ("submission", "release", "verified"))
    held, initial, final = sub["held_fields"], release["post_release_fields"], verified["scheduler_observations"][0]["fields"]
    for fields in (held, initial, final):
        assert fields["JobId"] == "22451" and fields["JobState"] == "PENDING"
        assert fields["Dependency"] == "afterany:22450(unfulfilled)"
        assert fields["NumCPUs"] == fields["CPUs/Task"] == "8"
        assert fields["MinMemoryNode"] == "64G" and fields["TimeLimit"] == "04:00:00"
        assert fields["Requeue"] == fields["Restarts"] == "0"
        assert fields["Comment"] == sub["request"]["sha256"] and fields["Command"] == sub["batch"]["path"]
    assert (held["Reason"], initial["Reason"], final["Reason"]) == ("JobHeldUser", "None", "Dependency")
    native, reviewer, conversion = [o["fields"] for o in verified["scheduler_observations"][1:]]
    assert native["JobId"] == "22444" and native["JobState"] == "RUNNING"
    assert reviewer["JobId"] == "22445" and reviewer["Dependency"] == "afterany:22444(unfulfilled)"
    assert conversion["JobId"] == "22450" and conversion["Dependency"] == "afterany:22445(unfulfilled)"
    assert all(o["command"][:3] == ["scontrol", "show", "job"] for o in verified["scheduler_observations"])


def test_frozen_scorer_preview_and_readonly_preflight_scope(receipts):
    preflight = receipts["preflight"]
    assert preflight["helper_sources_checked"] == 920
    assert preflight["scorer_runtime_records_checked"] == 704
    assert preflight["scorer_runtime_summed_record_bytes"] == 4746361652
    assert preflight["assessment_helper_records_checked"] == 12
    assert preflight["fas_protocol"]["newly_computed_pair_cap"] == 9000
    assert "unseeded" in preflight["fas_protocol"]["sampling"]
    command = preflight["command_preview"]
    assert "-resume" not in command
    assert command[command.index("--event_year") + 1] == "2020"
    assert command[command.index("--challenges_ids") + 1] == "GO EC VGNC SwissTrees TreeFam-A FAS"
    assert command[command.index("--participant_id") + 1] == "ohmm_qfo_full_native_p0_c1_r0"
    assert command[command.index("-work-dir") + 1] == preflight["output_namespaces"]["work"]
    runtime = preflight["runtime"]
    assert runtime["invocation"] != runtime["binary"]["path"]
    assert runtime["packages"] == dict(Bio="1.87", numpy="2.2.6", psutil="7.2.2")
    assert "unknown and potentially tool-dependent" in " ".join(preflight["limitations"])
