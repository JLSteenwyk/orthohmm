from copy import deepcopy

import pytest

from benchmark_tools.run_qfo_factorial_assessment import cell_label, verify_conversion_identity, verify_reuse_binding


def conversion(index):
    array = 21675 if index % 2 else 21674
    label = cell_label(index)
    accounting = f"JobID|JobIDRaw|State|ExitCode\n{array}_{index}|12345|COMPLETED|0:0\n"
    report = {"status": "cell_pairs_prepared_unscored", "accuracy_evaluated": False, "index": index,
              "cell": label, "participant": f"ohmm_qfo_factorial_{label}", "job_id": "12345",
              "array_job_id": str(array), "array_task_id": str(index), "retained_pairs": 8,
              "total_pairs": 10, "removed_mapping_pairs": 2,
              "semantics": "native phylogenetically inferred pairs" if index % 2 else "cross-species group-derived clique pairs"}
    return report, accounting


@pytest.mark.parametrize("index", range(8))
def test_all_cell_identities(index):
    report, accounting = conversion(index)
    assert verify_conversion_identity(report, index, accounting)["JobIDRaw"] == "12345"


@pytest.mark.parametrize("key,value", [("semantics", "RootHOG clique"), ("job_id", "wrong"),
    ("status", "running"), ("accuracy_evaluated", True), ("participant", "other"),
    ("removed_mapping_pairs", 5), ("retained_pairs", 11)])
def test_conversion_drift_rejected(key, value):
    report, accounting = conversion(3)
    report[key] = value
    with pytest.raises(ValueError):
        verify_conversion_identity(report, 3, accounting)


def test_live_failed_and_duplicate_jobs_rejected():
    report, accounting = conversion(2)
    for text in (accounting.replace("COMPLETED", "RUNNING"), accounting.replace("0:0", "1:0"),
                 accounting + accounting.splitlines()[1] + "\n"):
        with pytest.raises(ValueError):
            verify_conversion_identity(report, 2, text)


@pytest.mark.parametrize("index", [-1, 8, True, .5])
def test_invalid_cell_rejected(index):
    with pytest.raises(ValueError):
        cell_label(index)


def reuse():
    part = {"path": "/old", "sha256": "part", "bytes": 12}
    old = {"stage": "multipass_refined", "pairs": {"path": "/pairs", "sha256": "a"},
           "filtered_pairs": {"path": "/filtered", "sha256": "b"}, "total_pairs": 10,
           "retained_pairs": 8, "removed_mapping_pairs": 2, "partition": part}
    current = {**old, "candidate_partition": {**part, "path": "/candidate"},
               "reused_conversion": {"stage": "multipass_refined"}}
    return current, {"status": "admitted", "index": 1, "conversion": deepcopy(old)}


def test_exact_baseline_binding():
    current, native = reuse()
    verify_reuse_binding(current, native, 0)
    with pytest.raises(ValueError):
        verify_reuse_binding(current, native, 2)


@pytest.mark.parametrize("key", ["pairs", "filtered_pairs", "total_pairs", "candidate_partition", "reused_conversion"])
def test_baseline_changes_rejected(key):
    current, native = reuse()
    if key in {"pairs", "filtered_pairs"}:
        current[key] = {"path": "/different", "sha256": "other"}
    elif key == "total_pairs":
        current[key] = 11
    elif key == "candidate_partition":
        current[key] = {"sha256": "other", "bytes": 12}
    else:
        current[key] = {"stage": "other"}
    with pytest.raises(ValueError):
        verify_reuse_binding(current, native, 0)
