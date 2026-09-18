import pytest

from benchmark_tools.run_qfo_sequence_graph import validate_completion


def fixture():
    stages = ["initial", "multipass"]
    worker = dict(status="sequence_checked_worker_returned", calls=[
        dict(stage=s, status="checked", exit_code=0) for s in stages])
    metrics = dict(parameters=dict(profile_expansion=False),
        stages=[dict(label=s) for s in ("multipass", "multipass_refined")],
        counts=dict(genes=984137, species=78))
    variant = dict(expected_clustering_calls=stages, expected_stages=["multipass", "multipass_refined"],
                   expected_genes=984137, expected_species=78)
    return worker, metrics, variant


def test_complete():
    validate_completion(*fixture())


@pytest.mark.parametrize("problem", ["missing", "extra", "failed", "exit", "order", "profiles", "genes", "species", "stages"])
def test_completion_rejects(problem):
    worker, metrics, variant = fixture()
    if problem == "missing":
        worker["calls"].pop()
    elif problem == "extra":
        worker["calls"].append(dict(stage="profile_base", status="checked", exit_code=0))
    elif problem == "failed":
        worker["calls"][0]["status"] = "failed"
    elif problem == "exit":
        worker["calls"][0]["exit_code"] = 1
    elif problem == "order":
        worker["calls"].reverse()
    elif problem == "profiles":
        metrics["parameters"]["profile_expansion"] = True
    elif problem in ("genes", "species"):
        metrics["counts"][problem] = 1
    else:
        metrics["stages"].pop()
    with pytest.raises(ValueError):
        validate_completion(worker, metrics, variant)
