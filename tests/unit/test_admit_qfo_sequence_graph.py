import sys

import pytest

from benchmark_tools import admit_qfo_sequence_graph as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture(root):
    executor = root / "executor"
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    for base, names in ((executor, ("run_qfo_sequence_graph.py", *module.HELPERS)),
                        (launcher, ("replay_high_sensitivity.py",))):
        for name in names:
            path = base / "benchmark_tools" / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("# synthetic source\n")
    source = record(executor / "benchmark_tools/run_qfo_sequence_graph.py")
    plan_record = dict(path=str(root / "plan.json"), sha256="plan")
    variant = dict(native_command=["native"], requested_memory_gib=64,
        expected_clustering_calls=["initial", "multipass"], expected_stages=["multipass", "multipass_refined"],
        expected_genes=984137, expected_species=78)
    plan = dict(variants=dict(all_hits=variant), cwd=str(launcher), environment_overrides=dict(
        PYTHONPATH=str(launcher), PYTHONHASHSEED="0", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1"))
    parent = dict(status="sequence_checked_graph_complete_pending_admission", exit_code=0,
        job_id="job", executor_commit=module.EXECUTOR, variant="all_hits", accuracy_evaluated=False,
        source=source, plan=plan_record, helpers=[record(executor / "benchmark_tools" / n) for n in module.HELPERS],
        started_epoch=1, finished_epoch=2, worker_command=[sys.executable,
            str(executor / "benchmark_tools/run_qfo_sequence_graph.py"), "--root", str(root),
            "--plan", plan_record["path"], "--plan-sha256", "plan", "--variant", "all_hits", "--replay-worker"])
    scientific = record(launcher / "benchmark_tools/replay_high_sensitivity.py")
    worker = dict(status="sequence_checked_worker_returned", source=source, plan=plan_record,
        variant="all_hits", accuracy_evaluated=False, replay_source=scientific,
        calls=[dict(stage=s, status="checked", exit_code=0) for s in ("initial", "multipass")])
    replay = dict(source=scientific, cwd=str(launcher), command=["native"],
        stages=[dict(label=s) for s in ("multipass", "multipass_refined")], counts=dict(genes=984137, species=78),
        parameters=dict(accuracy_profile="high_sensitivity", cpm_resolution=.1, profile_expansion=False,
            profile_iterations=1, jackknife_profile_thresholds=False, jackknife_single_copy_profiles=False,
            profile_min_species=1, matrix="BLOSUM62", leiden_seed=4))
    scheduler = dict(JobIDRaw="job", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="32", ReqMem="64G")
    return parent, worker, replay, plan, plan_record, "all_hits", scheduler, root, executor


@pytest.mark.parametrize("problem", [None, "state", "memory", "job", "variant", "revision", "source", "helpers",
    "plan", "command", "cwd", "environment", "parameters", "score", "nan", "time_order", "worker_failed"])
def test_parent_gate(tmp_path, problem):
    args = fixture(tmp_path)
    parent, worker, replay, plan, _, _, scheduler, _, _ = args
    if problem == "state":
        scheduler["State"] = "RUNNING"
    elif problem == "memory":
        scheduler["ReqMem"] = "32G"
    elif problem == "job":
        parent["job_id"] = "other"
    elif problem == "variant":
        worker["variant"] = "top100"
    elif problem == "revision":
        parent["executor_commit"] = "other"
    elif problem == "source":
        worker["source"] = {}
    elif problem == "helpers":
        parent["helpers"].pop()
    elif problem == "plan":
        worker["plan"] = {}
    elif problem == "command":
        parent["worker_command"].pop()
    elif problem == "cwd":
        replay["cwd"] = "other"
    elif problem == "environment":
        plan["environment_overrides"]["OMP_NUM_THREADS"] = "32"
    elif problem == "parameters":
        replay["parameters"]["leiden_seed"] = 1
    elif problem == "score":
        replay["stages"][0]["official_orthobench"] = {}
    elif problem == "nan":
        parent["finished_epoch"] = float("nan")
    elif problem == "time_order":
        parent["finished_epoch"] = 0
    elif problem == "worker_failed":
        worker["calls"][0]["status"] = "failed"
    if problem:
        with pytest.raises(ValueError):
            module.validate_parent(*args)
    else:
        module.validate_parent(*args)


@pytest.mark.parametrize("problem", [None, "missing", "duplicate", "foreign", "count", "hash", "path"])
def test_partition_content(tmp_path, problem):
    directory = tmp_path / "replay"
    directory.mkdir()
    stages = []
    for label in ("multipass", "multipass_refined"):
        path = directory / ("orthogroups_" + label + ".txt")
        path.write_text("a b\nc\n")
        stages.append(dict(label=label, clusters=2, output=record(path)))
    path = directory / "orthogroups_multipass_refined.txt"
    if problem in ("missing", "duplicate", "foreign", "hash"):
        path.write_text({"missing": "a b\n", "duplicate": "a b\nc a\n", "foreign": "a b\nx\n", "hash": "a c\nb\n"}[problem])
        if problem != "hash":
            stages[1]["output"] = record(path)
    elif problem == "count":
        stages[1]["clusters"] = 3
    elif problem == "path":
        stages[1]["output"] = stages[0]["output"]
    if problem:
        with pytest.raises(ValueError):
            module.partitions(tmp_path, dict(stages=stages), {"a", "b", "c"})
    else:
        assert len(module.partitions(tmp_path, dict(stages=stages), {"a", "b", "c"})) == 2


def test_live_job_cannot_create_admission(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "sequence_evidence", lambda *args: (
        dict(variants=dict(all_hits=dict(requested_memory_gib=64))), {}, {}, {}))
    def fail(*args):
        raise ValueError("job not completed")
    monkeypatch.setattr(module, "completed", fail)
    destination = tmp_path / "admission.json"
    with pytest.raises(ValueError, match="not completed"):
        module.admit(tmp_path, tmp_path / "plan", "hash", "all_hits", "1", "hash", destination)
    assert not destination.exists()
