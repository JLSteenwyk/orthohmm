import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import run_qfo_cpm_phylogeny as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.prepare_orthobench_factorial import plan_cells
from tests.unit.test_run_qfo_parameter_phylogeny import allocation


def fixture(root, index=0):
    context = {"output_root": str(root / module.ARMS[index])}
    report = {"status": "cpm_candidates_admitted_unscored", "index": index, "arm": module.ARMS[index],
        "context": context, "accuracy_evaluated": False, "publication_ready": False,
        "candidate_arm": {"candidate_expansion": True,
            "seed_partition": {"path": str(Path(context["output_root"]) / "replay/orthogroups_profiles_refined.txt")},
            "candidate_partition": {"path": str(root / "candidates.txt")}, "membership_constraints": {"path": str(root / "merges.json")}}}
    return report, context


@pytest.mark.parametrize("index", [0, 1])
@pytest.mark.parametrize("problem", [None, "status", "index", "arm", "context", "accuracy", "ready", "seed", "expansion"])
def test_select_arm(tmp_path, index, problem):
    report, context = fixture(tmp_path, index)
    mutations = {"status": (report, "status", "failed"), "index": (report, "index", bool(index)),
        "arm": (report, "arm", "control"), "context": (report, "context", {}),
        "accuracy": (report, "accuracy_evaluated", True), "ready": (report, "publication_ready", True),
        "seed": (report["candidate_arm"], "seed_partition", {"path": "baseline.txt"}),
        "expansion": (report["candidate_arm"], "candidate_expansion", False)}
    if problem:
        target, key, value = mutations[problem]
        target[key] = value
        with pytest.raises(ValueError):
            module.select_arm(report, index, context)
    else:
        before = copy.deepcopy(report)
        arm = module.select_arm(report, index, context)
        assert arm["seed_partition"] == report["candidate_arm"]["seed_partition"]
        assert arm["label"] == module.ARMS[index] and report == before


@pytest.mark.parametrize("state", ["PENDING", "RUNNING", "FAILED"])
def test_live_scheduler_precedes_output_reads(tmp_path, monkeypatch, state):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        f"JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n22086_0|22087|{state}|0:0|00:01:00|bizon|2\n")
    with pytest.raises(ValueError, match="completed"):
        module.verify_sources(tmp_path, 0)


def test_allocation_precedes_sources(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="scheduled"):
        module.run(tmp_path, 0)


@pytest.mark.parametrize("outcome", ["success", "fresh", "failed", "exception", "changed"])
def test_execution_preserves_evidence_and_infers_own_tree(tmp_path, monkeypatch, outcome):
    allocation(monkeypatch)
    launcher = tmp_path / "launcher"
    launcher.mkdir()
    original = next(c for c in plan_cells(tmp_path / "baseline", tmp_path / "fastas",
        tmp_path / "executor/benchmark_tools/replay_phylogeny.py", 32) if c["label"] == "p1_c1_r1")
    report, context = fixture(tmp_path)
    arm = module.select_arm(report, 0, context)
    before = copy.deepcopy(original)
    verified = {"launcher": str(launcher), "prepared": str(tmp_path / "executor"), "original": original,
        "arm": arm, "environment": {}, "manifest": {"input_fastas": [], "environment_overrides": {"OMP_NUM_THREADS": "1"}},
        "checked_records": [], "admission_executor": str(tmp_path / "validator"), "admission": report}
    checks, launches = [], []
    def verify(*args):
        checks.append(args)
        return {} if outcome == "changed" and len(checks) > 1 else verified
    monkeypatch.setattr(module, "verify_sources", verify)
    monkeypatch.setattr(module, "native_command", lambda cell, *a: (cell["argv"], []))
    monkeypatch.setattr(module, "execution_environment", lambda *a: ({}, {}))
    def validate(command, **kwargs):
        Path(command[-1]).write_text(json.dumps({} if outcome == "fresh" else report))
    monkeypatch.setattr(module.subprocess, "run", validate)
    def execute(dataset, order, env, output, inputs, provenance):
        launches.append(dataset)
        assert Path.cwd() == launcher and order == ["candidate_cpm_low"]
        argv = dataset["methods"][order[0]]["argv"]
        assert argv[argv.index("--candidate-clusters") + 1] == arm["partition"]["path"]
        assert argv[argv.index("--membership-constraints") + 1] == arm["constraints"]["path"]
        assert argv[argv.index("--species-tree-mode") + 1] == "infer" and "--species-tree" not in argv
        assert env == {"OMP_NUM_THREADS": "1", "PYTHONPATH": str(launcher)}
        if outcome == "exception":
            raise RuntimeError("native interruption")
        return {"failed_methods": order if outcome == "failed" else []}
    monkeypatch.setattr(module, "execute", execute)
    cwd = Path.cwd()
    if outcome == "success":
        result = module.run(tmp_path, 0)
        assert result["status"] == "complete_pending_native_validation"
    else:
        with pytest.raises((ValueError, RuntimeError)):
            module.run(tmp_path, 0)
    destination = tmp_path / "benchmarks/results/qfo_cpm_phylogeny_v1/cpm_low"
    result = json.loads((destination / "postflight.json").read_text())
    assert result["status"] == ("complete_pending_native_validation" if outcome == "success" else "failed")
    assert result["accuracy_evaluated"] is False and result["native_outputs_validated"] is False
    assert len(launches) == (0 if outcome == "fresh" else 1)
    assert Path.cwd() == cwd and original == before
    with pytest.raises(FileExistsError):
        module.run(tmp_path, 0)


@pytest.mark.parametrize("problem", [None, "source", "revision", "file"])
def test_source_binding(tmp_path, monkeypatch, problem):
    executor = tmp_path / "benchmarks/work/publication_qfo_cpm_candidates_admission_v3"
    source = executor / "benchmark_tools/admit_qfo_cpm_candidates.py"
    source.parent.mkdir(parents=True)
    source.write_text("fixture")
    monkeypatch.setattr(module, "ADMISSION_SHA", "wrong" if problem == "source" else record(source)["sha256"])
    report, context = fixture(tmp_path)
    for item in (report["candidate_arm"][k] for k in ("candidate_partition", "membership_constraints", "seed_partition")):
        path = Path(item["path"])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture")
        item.update(record(path))
    report.update(source=record(source), checked_records=[])
    path = tmp_path / "benchmarks/work/qfo_cpm_candidates_admission_22086_0.json"
    path.write_text(json.dumps(report))
    if problem == "file":
        Path(report["candidate_arm"]["seed_partition"]["path"]).write_text("changed")
    monkeypatch.setattr(module, "corrected_evidence", lambda *a: ({}, {}, {}, {}))
    context["checked_records"] = []
    # Persist the exact context after adding its source records.
    report["context"] = context
    path.write_text(json.dumps(report))
    monkeypatch.setattr(module, "evidence", lambda *a: context)
    monkeypatch.setattr(module, "verify_baseline", lambda *a: ({}, {}, tmp_path, tmp_path, {}, []))
    def check_output(command, **kwargs):
        if command[0] == "sacct":
            return "JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n22086_0|22087|COMPLETED|0:0|00:01:00|bizon|2\n"
        return "wrong" if problem == "revision" else module.ADMISSION_COMMIT
    monkeypatch.setattr(module.subprocess, "check_output", check_output)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    if problem:
        with pytest.raises(ValueError):
            module.verify_sources(tmp_path, 0)
    else:
        result = module.verify_sources(tmp_path, 0)
        assert result["admission"] == report
        assert result["arm"]["label"] == "cpm_low"
