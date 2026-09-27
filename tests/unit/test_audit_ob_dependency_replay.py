import json

import numpy as np
import pytest

from benchmark_tools.audit_ob_dependency_replay import (
    ENV, candidate_readback, compare_graphs, validate_execution,
)
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def execution_fixture():
    arm = dict(python="python", directory="arm")
    plan = dict(source={"path": "driver"}, repo="repo", checkpoint="checkpoint", inputs="inputs", cpu=32)
    execution = dict(returncode=0, timed_out=False, attempts=1, job_id="22320", env=dict(ENV),
        command=["python", "-I", "driver", "--driver", "arm", "--repo", "repo", "--checkpoint", "checkpoint",
                 "--inputs", "inputs", "--cpu", "32"])
    return execution, arm, plan


def test_execution_matches():
    validate_execution(*execution_fixture())


@pytest.mark.parametrize("key,value", [("returncode", 1), ("timed_out", True), ("attempts", 2),
    ("job_id", "another"), ("env", {}), ("command", [])])
def test_execution_rejects_changed_contract(key, value):
    execution, arm, plan = execution_fixture()
    execution[key] = value
    with pytest.raises(ValueError):
        validate_execution(execution, arm, plan)


@pytest.mark.parametrize("change", [None, "weights", "endpoints", "invalid", "extra"])
def test_graph_comparison(tmp_path, change):
    values = dict(sources=np.array([0, 0]), targets=np.array([1, 2]), weights=np.array([1., 2.]))
    left, right = tmp_path / "a.npz", tmp_path / "b.npz"
    np.savez(left, **values)
    if change == "weights":
        values["weights"][1] = 3.
    elif change == "endpoints":
        values["sources"][1] = 1
    elif change == "invalid":
        values["weights"][1] = np.nan
    elif change == "extra":
        values["extra"] = np.array([0])
    np.savez(right, **values)
    if change in {"invalid", "extra"}:
        with pytest.raises(ValueError):
            compare_graphs(left, right, list("abc"))
    else:
        result = compare_graphs(left, right, list("abc"))
        assert result["arrays_equal"] == (change is None)
        assert result["endpoints_equal"] == (change != "endpoints")


def candidate_fixture(tmp_path):
    work = tmp_path / "replay/orthohmm_working_res"
    work.mkdir(parents=True)
    seed = tmp_path / "seeds.txt"
    seed.write_text("a b\nc\nd\n")
    for name in ("orthohmm_edges_clustered.txt", "phylogeny_candidate_superfamilies.txt"):
        (work / name).write_text("a b c\nd\n")
    (work / "phylogeny_candidate_seeds.tsv").write_text(
        "candidate_family\tseed_families\nFamily0000000\tSeed0000000,Seed0000001\nFamily0000001\tSeed0000002\n")
    (work / "phylogeny_candidate_merges.json").write_text(json.dumps([
        dict(source_genes=["c"], target_genes=["a", "b"], source_size=1, target_size=2,
             source_seed_families=1, target_seed_families=1)]))
    stage = dict(candidates=dict(profile="satellite_v2", membership_policy="high_confidence_pair",
        candidate_checkpoint=str(work / "phylogeny_candidate_superfamilies.txt"),
        seed_sidecar=str(work / "phylogeny_candidate_seeds.tsv"),
        merge_trace_sidecar=str(work / "phylogeny_candidate_merges.json"),
        seed_families=3, candidate_families=2, merges=1))
    return stage, record(seed), work


def test_candidate_adapter_preserves_original_evidence(tmp_path):
    stage, seed, work = candidate_fixture(tmp_path)
    before = [record(p) for p in sorted(work.iterdir())]
    result = candidate_readback(tmp_path, stage, seed, set("abcd"))
    assert result["candidate_families"] == 2 and result["merges"] == 1
    assert sorted(result["original_records"], key=lambda r: r["path"]) == before
    assert before == [record(p) for p in sorted(work.iterdir())]
    assert "checked_records" not in result


@pytest.mark.parametrize("problem", ["path", "sidecar", "trace", "split", "count"])
def test_candidate_adapter_rejects_inconsistent_content(tmp_path, problem):
    stage, seed, work = candidate_fixture(tmp_path)
    if problem == "path":
        stage["candidates"]["seed_sidecar"] = "wrong"
    elif problem == "sidecar":
        (work / "phylogeny_candidate_seeds.tsv").write_text("wrong\n")
    elif problem == "trace":
        (work / "phylogeny_candidate_merges.json").write_text("[]")
    elif problem == "split":
        (work / "orthohmm_edges_clustered.txt").write_text("a c\nb d\n")
    else:
        stage["candidates"]["merges"] = 2
    with pytest.raises(ValueError):
        candidate_readback(tmp_path, stage, seed, set("abcd"))
