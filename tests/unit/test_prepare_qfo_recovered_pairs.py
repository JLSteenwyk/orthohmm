import pytest

from benchmark_tools.prepare_qfo_recovered_pairs import stage_outputs, STAGES, prepare


def fixture():
    stages = [{"label": s, "output": {"path": s}} for s in STAGES]
    return {"status": "checked_full_replay_recovered_verified", "accuracy_evaluated": False,
            "run_version": "v2", "recovery": {"inference_rerun": False}, "native_replay": {"stages": stages},
            "coverage": [{"stage": s["label"], "observed": s["output"].copy(), "partition_equal": True} for s in stages]}


@pytest.mark.parametrize("problem", [None, "status", "scored", "version", "rerun", "missing", "order", "coverage", "partition"])
def test_all_four_admitted_stages_required(problem):
    report = fixture()
    if problem == "status":
        report["status"] = "failed"
    elif problem == "scored":
        report["accuracy_evaluated"] = True
    elif problem == "version":
        report["run_version"] = "v1"
    elif problem == "rerun":
        report["recovery"]["inference_rerun"] = True
    elif problem == "missing":
        report["native_replay"]["stages"].pop()
    elif problem == "order":
        report["native_replay"]["stages"].reverse()
    elif problem == "coverage":
        report["coverage"].pop()
    elif problem == "partition":
        report["coverage"][0]["observed"] = {"path": "other"}
    if problem:
        with pytest.raises(ValueError):
            stage_outputs(report)
    else:
        assert len(stage_outputs(report)) == 4


def test_no_overwrite(tmp_path):
    with pytest.raises(FileExistsError):
        prepare(tmp_path, tmp_path)
