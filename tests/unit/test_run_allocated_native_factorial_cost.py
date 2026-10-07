from copy import deepcopy
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import run_allocated_native_factorial_cost as controller
from benchmark_tools.native_factorial_adapter import factors
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_allocated_threadripper_scaling import placement_value, ALLOWED, CGROUP


def ready():
    value = placement_value()
    return dict(pid=1, cgroup=CGROUP, placement=value["bound"], allocated_placement=value)


def test_actual_native_child_inherits_selected_cores():
    r = ready()
    current = dict(deepcopy(r["placement"]), pid=2)
    assert controller.native_placement(r, current, 42, 1) == ALLOWED
    assert current["affinity"] != list(range(32))


@pytest.mark.parametrize("change", ["same_pid", "parent", "cgroup", "cpus", "job", "ram", "ancestors", "ready"])
def test_native_inheritance_refusals(change):
    r = ready()
    current = dict(deepcopy(r["placement"]), pid=2)
    parent = 1
    if change == "same_pid": current["pid"] = 1
    elif change == "parent": parent = 999
    elif change == "cgroup": current["cgroup"] = CGROUP.replace("step_0", "step_1")
    elif change == "cpus": current["affinity"] = list(range(32))
    elif change == "job": current["slurm"]["SLURM_JOB_ID"] = "43"
    elif change == "ram":
        for row in current["ancestors"]: row["memory.max"] = str(64*1024**3)
    elif change == "ancestors": current["ancestors"][0]["cpu.max"] = "100000 100000"
    else: r["placement"]["pid"] = 999
    with pytest.raises(ValueError):
        controller.native_placement(r, current, 42, parent)


@pytest.fixture
def native_case(tmp_path, monkeypatch):
    root = tmp_path / "run_10"
    (root / "measurement").mkdir(parents=True)
    (root / "input").mkdir()
    (root / "measurement/ready.json").write_text(json.dumps(ready()))
    amendment_path = tmp_path / "amendment.json"
    amendment_path.write_text("{}")
    amendment_ref = record(amendment_path)
    baseline_path = tmp_path / "baseline.json"
    core = tmp_path / "core"
    (core / "orthohmm").mkdir(parents=True)
    source = core / "orthohmm/orthohmm.py"
    source.write_text("# Synthetic entrypoint stub, not the frozen scientific pipeline.\n")
    baseline = dict(core_root=str(core), environment_overrides={"PYTHONHASHSEED": "0"},
        tool_entrypoints=dict(orthohmm_python=dict(absolute_path=str(Path(sys.executable).absolute())),
            mafft=dict(absolute_path="/frozen/mafft"), FastTree=dict(absolute_path="/frozen/FastTree")))
    baseline_path.write_text(json.dumps(baseline))
    run = dict(index=10, cell="p1_c0_r1", output_root=str(root), native_order=["s0.fa"], genes=8)
    plan = dict(baseline=record(baseline_path), runs=[{}]*10+[run])
    execution = dict(allowed_indices=[10,11,12], historical_plan=dict(path="/plan", bytes=1, sha256="a"*64))
    monkeypatch.setattr(controller, "amendment", lambda ref: (execution, plan))
    monkeypatch.setenv("SLURM_JOB_ID", "42")
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    monkeypatch.setattr(controller.sys, "dont_write_bytecode", True)
    monkeypatch.setattr(controller.sys, "path", list(sys.path))
    monkeypatch.setattr(controller.sys, "pycache_prefix", str(tmp_path / "absent_cache"))
    monkeypatch.setattr(controller.os, "getppid", lambda: 1)
    from benchmark_tools import probe_threadripper_allocation
    current = dict(deepcopy(ready()["placement"]), pid=2)
    monkeypatch.setattr(probe_threadripper_allocation, "inspect", lambda: current)
    module = SimpleNamespace(__file__=str(source), fetch_fasta_files=lambda path: run["native_order"],
        StopStep=SimpleNamespace(infer="infer"), SubstitutionMatrix=SimpleNamespace(blosum62="BLOSUM62"))
    original_import = controller.importlib.import_module
    monkeypatch.setattr(controller.importlib, "import_module", lambda name, *a, **k:
        module if name == "orthohmm.orthohmm" else original_import(name, *a, **k))
    calls = []
    state = dict(error=None, bad_metrics=False)
    flags = factors(run["cell"])
    def adapted(**kwargs):
        calls.append(kwargs)
        if state["error"]:
            raise state["error"]
        stages = {"search", "edge_thresholds", "network_edges", "clustering", "refinement",
            "orthogroup_materialization", "profile_expansion", "phylogeny"}
        metrics = dict(status="complete", metadata=dict(native_factorial=flags), stages=dict.fromkeys(stages,{}),
            counts=dict(genes=7 if state["bad_metrics"] else 8, phylogeny_checkpoint_hits=0,
                        phylogeny_species_tree_checkpoint_hit=False))
        (root / "metrics.json").write_text(json.dumps(metrics))
    monkeypatch.setattr(controller, "entrypoint", lambda *a: (adapted, flags))
    yield dict(root=root, amendment=amendment_ref, calls=calls, state=state, run=run, flags=flags,
               execution=execution, plan=plan, current=current)


def test_native_receipt_and_frozen_arguments(native_case):
    d = native_case
    controller.native(d["amendment"], 10)
    receipt = json.loads((d["root"] / "native_execution.json").read_text())
    assert receipt["schema"] == "allocated_native_factorial_execution_v1"
    assert receipt["status"] == "native_factorial_completed_pending_output_review"
    assert receipt["native_cpu_ids"] == ALLOWED and receipt["parent_pid"] == 1
    assert receipt["amendment"] == d["amendment"]
    assert receipt["automatic_retry"] is False and receipt["accuracy_evaluated"] is False
    args = d["calls"][0]
    assert args["cpu"] == 32 and args["threads_per_worker"] == 4
    assert args["accuracy_profile"] == "high_sensitivity" and args["search_mode"] == "builtin"
    assert args["phylogeny"] == "reconcile" and args["phylogeny_candidates"] == "seed"
    assert args["clustering"] == "leiden" and args["cpm_resolution"] == .1
    assert args["evalue_threshold"] == .0001 and args["substitution_matrix"] == "BLOSUM62"


@pytest.mark.parametrize("index", [9, 13, True])
def test_completed_or_invalid_identity_never_executes(native_case, index):
    with pytest.raises(ValueError):
        controller.native(native_case["amendment"], index)
    assert native_case["calls"] == []


@pytest.mark.parametrize("bad", ["kernel", "metrics", "reused_output", "placement", "cache"])
def test_native_failures_retained_no_retry(native_case, bad):
    d = native_case
    if bad == "kernel": d["state"]["error"] = ValueError("synthetic kernel failure")
    elif bad == "metrics": d["state"]["bad_metrics"] = True
    elif bad == "reused_output": (d["root"] / "native").mkdir()
    elif bad == "placement": d["current"]["affinity"] = list(range(32))
    else: Path(sys.pycache_prefix).mkdir()
    with pytest.raises((ValueError, FileExistsError)):
        controller.native(d["amendment"], 10)
    path = d["root"] / "native_execution.json"
    if bad in {"kernel", "metrics"}:
        result = json.loads(path.read_text())
        assert result["status"] == "native_factorial_failed" and result["automatic_retry"] is False
    else:
        assert not path.exists() and d["calls"] == []
