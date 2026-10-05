import copy
import hashlib
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import export_native_factorial_progress as report


def write(path, value):
    data = (json.dumps(value, sort_keys=True) + "\n").encode()
    path.write_bytes(data)
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


@pytest.fixture
def sample(tmp_path):
    cells = ["p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p0_c1_r1", "p1_c0_r1", "p1_c1_r0"]
    runs = [{"index": i, "dataset": "orthobench" if i < 6 else "qfo_corrected",
             "cell": (cells + cells + ["p1_c1_r1"])[i], "repeat": 0} for i in range(13)]
    plan = {"runs": runs}
    plan_ref = write(tmp_path / "plan.json", plan)
    review = {"schema": "native_factorial_terminal_review_v1", **runs[1], "job_id": 22428,
              "terminal_reviewed": True, "primary_resources_replayed": True,
              "shared_host_resources_reviewed": True, "native_outputs_validated": True,
              "execution_scope": "shared_host_matched_resources", "uncontended_timing": False,
              "resource_scopes": report.SCOPES, "plan": plan_ref,
              "status": "native_success", "scheduler_state": "COMPLETED", "scheduler_exit_code": "0:0",
              "resources": {"wall_seconds": 100.5, "cpu_seconds": 3000.1, "peak_memory_bytes": 2**30},
              "whole_run_maximum_foreign_average_cores": 43.1}
    review_ref = write(tmp_path / "review.json", review)
    score = {"schema": "native_factorial_orthobench_score_v1", "status": "terminal_native_orthobench_scored",
             "dataset": "OrthoBench", "index": 1, "cell": cells[1], "job_id": 22428,
             "terminal_review": review_ref, "native_outputs_validated": True, "accuracy_evaluated": True,
             "plan": plan_ref, "resource_observation": review["resources"], "resource_scopes": report.SCOPES,
             "canonical_partition_comparison": {"partition_equal": False}, "score_percent": {"refogs": 70},
             "reference_families": 70, "reference_genes": 1944, "covered_reference_genes": 1900,
             "development_exposed": True, "independent_validation": False,
             "score_fraction": {"f_score": .7, "precision": .8, "recall": .6}}
    score_ref = write(tmp_path / "score.json", score)
    return tmp_path, plan, plan_ref, review, review_ref, score, score_ref


def collect(sample, with_score=True):
    _, _, plan_ref, _, review_ref, _, score_ref = sample
    args = [review_ref["path"], review_ref["sha256"],
            score_ref["path"] if with_score else "-", score_ref["sha256"] if with_score else "-"]
    return report.collect(plan_ref["path"], plan_ref["sha256"], [args])


def change(sample, target, key, value):
    offset = {"plan": 1, "review": 3, "score": 5}[target]
    sample[offset][key] = value
    sample[offset + 1].update(write(Path(sample[offset + 1]["path"]), sample[offset]))


def test_success_keeps_unavailable_distinct_from_zero(sample):
    result = collect(sample)
    assert len(result["rows"]) == 13
    assert result["rows"][0]["wall_seconds"] is None
    row = result["rows"][1]
    assert row["outcome"] == "native_success" and row["f1"] == .7
    assert row["reference_gene_coverage"] == 1900 / 1944
    assert row["partition_equal"] is False
    assert row["maximum_foreign_average_cores"] == 43.1
    text = report.markdown(result)
    assert report.DISCLOSURE in text and "70.0000" in text and "Unavailable" in text
    assert result["new_scoring_or_admission"] is False
    assert result["publication_ready"] is False


def test_unscored_success_does_not_inherit_accuracy(sample):
    result = collect(sample, False)
    assert result["rows"][1]["f1"] is None
    assert result["rows"][1]["wall_seconds"] == 100.5


@pytest.mark.parametrize("key,value", [
    ("terminal_reviewed", False), ("primary_resources_replayed", False),
    ("shared_host_resources_reviewed", False), ("uncontended_timing", True),
    ("native_outputs_validated", False), ("scheduler_state", "RUNNING"),
    ("scheduler_exit_code", "1:0"), ("status", "queued"),
    ("index", True), ("index", -1), ("cell", "p1_c0_r0"),
    ("resource_scopes", {}), ("resources", None),
    ("whole_run_maximum_foreign_average_cores", None),
    ("whole_run_maximum_foreign_average_cores", float("nan")),
])
def test_bad_review_refused(sample, key, value):
    change(sample, "review", key, value)
    with pytest.raises(ValueError):
        collect(sample, False)


@pytest.mark.parametrize("key,value", [
    ("native_outputs_validated", False), ("accuracy_evaluated", False),
    ("independent_validation", True), ("development_exposed", False),
    ("reference_families", 69), ("covered_reference_genes", 2000),
    ("resource_observation", {}), ("job_id", 123),
    ("score_fraction", {"f_score": float("inf"), "precision": .8, "recall": .6}),
    ("score_fraction", {"f_score": -.1, "precision": .8, "recall": .6}),
])
def test_bad_score_refused(sample, key, value):
    change(sample, "score", key, value)
    with pytest.raises(ValueError):
        collect(sample)


def test_duplicate_and_checksum_refused(sample):
    _, _, plan_ref, _, review_ref, _, _ = sample
    attempt = [review_ref["path"], review_ref["sha256"], "-", "-"]
    with pytest.raises(ValueError, match="duplicate"):
        report.collect(plan_ref["path"], plan_ref["sha256"], [attempt, attempt])
    with pytest.raises(ValueError, match="checksum"):
        report.collect(plan_ref["path"], "0" * 64, [])


def test_recovered_failure_never_becomes_success(sample):
    tmp, plan, plan_ref, review, _, _, _ = sample
    old_plan_ref = copy.deepcopy(plan_ref)
    review.update({**plan["runs"][0], "job_id": 22427, "status": "native_failure_retained",
                   "native_outputs_validated": False, "scheduler_state": "FAILED", "scheduler_exit_code": "1:0"})
    review_ref = write(tmp / "failed_review.json", review)
    data = {"schema": "recovered_native_orthobench_score_v1", "cell": review["cell"],
            "dataset": "OrthoBench", "timing_success_established": False,
            "original_partition_comparison": {"partition_equal": True}, "score": {"refogs": 70},
            "development_exposed": True, "independent_validation": False,
            "score_fraction": {"f_score": .69, "precision": .78, "recall": .62}}
    data_ref = write(tmp / "recovered_score.json", data)
    recovery = {"schema": "native_factorial_receipt_failure_recovery_v1",
                "status": "scientific_outputs_recovered_from_failed_wrapper", "index": 0,
                "cell": review["cell"], "job_id": 22427, "terminal_review": review_ref,
                "native_outputs_validated": True, "accuracy_evaluated": True,
                "scheduler_success": False, "native_command_success": False,
                "timing_success_established": False, "resources": review["resources"],
                "resource_scopes": report.SCOPES, "score": data_ref}
    recovery_ref = write(tmp / "recovery.json", recovery)
    plan["history_adoption"] = {"prior_plan": old_plan_ref, "prior_review": review_ref,
                                "recovery": recovery_ref}
    plan_ref = write(tmp / "new_plan.json", plan)
    result = report.collect(plan_ref["path"], plan_ref["sha256"],
                            [[review_ref["path"], review_ref["sha256"], recovery_ref["path"], recovery_ref["sha256"]]])
    row = result["rows"][0]
    assert row["outcome"] == "failed_wrapper_science_recovered" and row["f1"] == .69
    assert row["reference_gene_coverage"] is None
    assert row["wall_seconds"] == 100.5
    assert "not clean-success timings" in report.markdown(result)


def test_pre_native_failure_has_no_imputed_resources(sample):
    review = sample[3]
    review.update(status="native_failure_retained", scheduler_state="FAILED", scheduler_exit_code="1:0",
                  native_outputs_validated=False, resources=None, whole_run_maximum_foreign_average_cores=None)
    change(sample, "review", "status", "native_failure_retained")
    row = collect(sample, False)["rows"][1]
    assert row["wall_seconds"] is None and row["maximum_foreign_average_cores"] is None


def test_qfo_refuses_orthobench_score(sample):
    change(sample, "review", "index", 6)
    change(sample, "review", "dataset", "qfo_corrected")
    change(sample, "review", "cell", "p0_c0_r0")
    with pytest.raises(ValueError, match="QfO"):
        collect(sample)


def test_report_file_symlink_refused(sample):
    tmp, _, ref, _, _, _, _ = sample
    link = tmp / "link.json"
    link.symlink_to(ref["path"])
    with pytest.raises(ValueError, match="Nonregular"):
        report.collect(link, ref["sha256"], [])


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1, True, "3000"])
def test_invalid_resources_refused(sample, value):
    resources = dict(sample[3]["resources"], cpu_seconds=value)
    change(sample, "review", "resources", resources)
    with pytest.raises(ValueError, match="numeric"):
        collect(sample, False)


def test_unknown_plan_history_not_silently_adopted(sample):
    change(sample, "review", "plan", {"bytes": 10, "sha256": "0" * 64})
    with pytest.raises(ValueError, match="outside plan"):
        collect(sample, False)


def test_qfo_resources_can_be_reported_without_synthetic_f1(sample):
    change(sample, "review", "index", 6)
    change(sample, "review", "dataset", "qfo_corrected")
    change(sample, "review", "cell", "p0_c0_r0")
    row = collect(sample, False)["rows"][6]
    assert row["wall_seconds"] == 100.5 and row["f1"] is None


def test_main_writes_all_rows_and_refuses_overwrite(sample, monkeypatch):
    tmp, _, plan_ref, _, review_ref, _, score_ref = sample
    output = tmp / "output"
    argv = ["export", "--plan", plan_ref["path"], "--plan-sha256", plan_ref["sha256"],
            "--attempt", review_ref["path"], review_ref["sha256"], score_ref["path"], score_ref["sha256"],
            "--output", str(output)]
    monkeypatch.setattr(sys, "argv", argv)
    report.main()
    result = json.loads((output / "report.json").read_text())
    assert len((output / "rows.tsv").read_text().splitlines()) == 14
    assert (output / "report.md").read_text() == report.markdown(result)
    assert result["source"]["sha256"] == hashlib.sha256(Path(report.__file__).read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="Output already exists"):
        report.main()


def test_incomplete_optional_score_arguments_refused(sample):
    _, _, plan_ref, _, review_ref, _, _ = sample
    with pytest.raises(ValueError, match="optional score"):
        report.collect(plan_ref["path"], plan_ref["sha256"],
                       [[review_ref["path"], review_ref["sha256"], "-", "0" * 64]])
