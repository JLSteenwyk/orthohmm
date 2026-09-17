import pytest

from benchmark_tools.admit_qfo_checked_full_replay import CALLS, STAGES, check_inventory, admit


def fixture():
    parent = {"status": "full_checked_replay_complete_unscored", "job_id": "21329", "exit_code": 0, "accuracy_evaluated": False}
    worker = {"status": "checked_full_replay_returned", "accuracy_evaluated": False, "calls": [
        {"index": i, "stage": stage, "status": "checked", "exit_code": 0, "accuracy_evaluated": False} for i, stage in enumerate(CALLS)]}
    replay = {"stages": [{"label": label} for label in STAGES], "counts": {"genes": 976504, "profiles_built": 12},
              "parameters": {"accuracy_profile": "high_sensitivity", "cpm_resolution": .1, "profile_expansion": True,
                "profile_iterations": 1, "jackknife_profile_thresholds": False,
                "jackknife_single_copy_profiles": False, "profile_min_species": 1, "matrix": "BLOSUM62", "leiden_seed": 4}}
    return parent, worker, replay


@pytest.mark.parametrize("problem", [None, "parent_running", "job", "worker_running", "missing_call", "call_order",
    "call_failed", "missing_stage", "stage_order", "missing_profiles", "zero_profiles", "genes", "parameters", "scored"])
def test_exact_unscored_inventory(problem):
    parent, worker, replay = fixture()
    if problem == "parent_running":
        parent["status"] = "running"
    elif problem == "job":
        parent["job_id"] = "21328"
    elif problem == "worker_running":
        worker["status"] = "running"
    elif problem == "missing_call":
        worker["calls"].pop()
    elif problem == "call_order":
        worker["calls"].reverse()
    elif problem == "call_failed":
        worker["calls"][0]["exit_code"] = 1
    elif problem == "missing_stage":
        replay["stages"].pop()
    elif problem == "stage_order":
        replay["stages"].reverse()
    elif problem == "missing_profiles":
        del replay["counts"]["profiles_built"]
    elif problem == "zero_profiles":
        replay["counts"]["profiles_built"] = 0
    elif problem == "genes":
        replay["counts"]["genes"] -= 1
    elif problem == "parameters":
        replay["parameters"]["leiden_seed"] = 9
    elif problem == "scored":
        replay["stages"][0]["official_orthobench"] = {}
    if problem:
        with pytest.raises(ValueError):
            check_inventory(parent, worker, replay)
    else:
        check_inventory(parent, worker, replay)


def test_refuse_reusing_admission(tmp_path):
    with pytest.raises(FileExistsError):
        admit(tmp_path, tmp_path, "unused")


@pytest.mark.parametrize("problem", [None, "old_job", "missing_environment", "wrong_inherited", "wrong_override"])
def test_v2_requires_pinned_job_and_child_thread_isolation(problem):
    parent, worker, replay = fixture()
    parent["job_id"] = "21333"
    for row in worker["calls"]:
        overrides = {k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
        row["thread_environment"] = {"child_overrides": overrides,
            "inherited": {**overrides, "OMP_NUM_THREADS": "32" if row["index"] >= 2 else "1"}}
    if problem == "old_job":
        parent["job_id"] = "21329"
    elif problem == "missing_environment":
        del worker["calls"][2]["thread_environment"]
    elif problem == "wrong_inherited":
        worker["calls"][2]["thread_environment"]["inherited"]["OMP_NUM_THREADS"] = "1"
    elif problem == "wrong_override":
        worker["calls"][2]["thread_environment"]["child_overrides"]["OMP_NUM_THREADS"] = "32"
    if problem:
        with pytest.raises(ValueError):
            check_inventory(parent, worker, replay, "v2")
    else:
        check_inventory(parent, worker, replay, "v2")


def test_unknown_run_version_rejected():
    with pytest.raises(ValueError, match="Unknown"):
        check_inventory(*fixture(), version="v3")
