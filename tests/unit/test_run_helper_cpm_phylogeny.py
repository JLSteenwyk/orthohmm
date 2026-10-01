import copy
import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import run_helper_cpm_phylogeny as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.prepare_orthobench_factorial import plan_cells


def report_fixture(root):
    directory = root / "benchmarks/results/qfo_cpm_helper_recovered_candidates_v1/candidate/orthohmm_working_res"
    arm = {"candidate_expansion": True,
        "candidate_partition": {"path": str(directory / "orthohmm_edges_clustered.txt")},
        "membership_constraints": {"path": str(directory / "phylogeny_candidate_merges.json")},
        "seed_partition": {"path": str(root / "benchmarks/results/qfo_cpm_checkpoint_recovery_v1/orthogroups_profiles_refined.txt")},
        "expansion": {"profile": "satellite_v2", "membership_policy": "high_confidence_pair"}}
    return {"status": "cpm_helper_recovered_candidates_admitted_unscored", "candidate_admitted": True,
        "arm": "cpm_high", "index": 1, "seed_handoff": "explicit_helper_runtime_seed_amendment",
        "accuracy_evaluated": False, "downstream_admitted": False, "publication_ready": False,
        "candidate_arm": arm, "verification": dict(genes=984137, seed_groups=390845,
            candidate_groups=346866, reconstructed_merges=43979)}


@pytest.mark.parametrize("problem", [None, "status", "admitted", "arm", "index", "handoff", "accuracy",
    "downstream", "ready", "counts", "expansion", "seed", "partition", "constraints", "profile", "policy"])
def test_exact_recovered_arm(tmp_path, problem):
    report = report_fixture(tmp_path)
    arm = report["candidate_arm"]
    mutations = {"status": (report, "status", "cpm_candidates_admitted_unscored"),
        "admitted": (report, "candidate_admitted", 1), "arm": (report, "arm", "cpm_low"),
        "index": (report, "index", True), "handoff": (report, "seed_handoff", "historical"),
        "accuracy": (report, "accuracy_evaluated", True), "downstream": (report, "downstream_admitted", True),
        "ready": (report, "publication_ready", True), "counts": (report, "verification", {}),
        "expansion": (arm, "candidate_expansion", 1), "seed": (arm, "seed_partition", {"path": "old_seed"}),
        "partition": (arm, "candidate_partition", {"path": "old_candidates"}),
        "constraints": (arm, "membership_constraints", {"path": "old_constraints"}),
        "profile": (arm["expansion"], "profile", "satellite"),
        "policy": (arm["expansion"], "membership_policy", "unconstrained")}
    if problem:
        target, key, value = mutations[problem]
        target[key] = value
        with pytest.raises(ValueError):
            module.select_arm(tmp_path, report)
    else:
        before = copy.deepcopy(report)
        selected = module.select_arm(tmp_path, report)
        assert selected == dict(label="cpm_high", partition=arm["candidate_partition"],
            constraints=arm["membership_constraints"], seed_partition=arm["seed_partition"])
        assert report == before


def accounting(state="COMPLETED", exit_code="0:0", cpus="2", memory="64G", node="bizon"):
    return ("JobID|JobIDRaw|State|ExitCode|Elapsed|AllocCPUS|ReqMem|NodeList\n"
            f"22386|22386|{state}|{exit_code}|00:01:14|{cpus}|{memory}|{node}\n")


@pytest.mark.parametrize("values", [{"state": "RUNNING"}, {"state": "FAILED"}, {"exit_code": "1:0"},
    {"cpus": "32"}, {"memory": "192G"}, {"node": "dgx"}])
def test_parent_allocation_precedes_file_reads(tmp_path, monkeypatch, values):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: accounting(**values))
    with pytest.raises(ValueError, match="completed"):
        module.verify_candidates(tmp_path)


@pytest.mark.parametrize("text", ["", accounting() + accounting().splitlines()[1] + "\n"])
def test_missing_or_duplicate_parent(text):
    with pytest.raises(ValueError):
        module.completed(text)


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return record(path)


def source_fixture(root, monkeypatch, problem):
    executor = root / module.EXECUTOR
    source = executor / "benchmark_tools/admit_helper_cpm_candidates.py"
    source.parent.mkdir(parents=True)
    source.write_text("frozen admission source\n")
    report = report_fixture(root)
    records = []
    for key in ("seed_partition", "candidate_partition", "membership_constraints"):
        path = Path(report["candidate_arm"][key]["path"])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("frozen content\n")
        report["candidate_arm"][key] = record(path)
        records.append(record(path))
    report.update(source=record(source), checked_records=[record(source), *records])
    admission = write_json(root / module.ADMISSION, report)
    extra = root / "submission.json"
    extra.write_text("{}\n")
    readback = {"status": "recovered_candidate_admission_independently_read_back_unscored",
        "admission": admission, "candidate_arm": report["candidate_arm"], "verification": report["verification"],
        "candidate_admitted": True, "accuracy_evaluated": False, "downstream_admitted": False,
        "publication_ready": False, "source_revision": module.COMMIT,
        "submission": record(extra), "log": record(extra), "time": record(extra),
        "source_git_bindings": [{**record(source), "git_revision": module.COMMIT}]}
    if problem == "readback_semantics":
        readback["downstream_admitted"] = True
    if problem == "git_revision":
        readback["source_git_bindings"][0]["git_revision"] = "wrong"
    if problem == "admission_binding":
        readback["admission"] = record(extra)
    if problem == "source_binding":
        report["source"] = record(extra)
    if problem == "conflict":
        report["checked_records"].append({**record(source), "sha256": "wrong"})
    if problem in ("conflict", "source_binding"):
        admission = write_json(root / module.ADMISSION, report)
        readback["admission"] = admission
    readback_record = write_json(root / module.READBACK, readback)
    monkeypatch.setattr(module, "ADMISSION_SHA", "wrong" if problem == "admission_pin" else admission["sha256"])
    monkeypatch.setattr(module, "READBACK_SHA", "wrong" if problem == "readback_pin" else readback_record["sha256"])
    monkeypatch.setattr(module, "SOURCE_SHA", "wrong" if problem == "source_pin" else record(source)["sha256"])
    def check_output(command, **kwargs):
        if command[0] == "sacct":
            return accounting()
        if command[-2] == "show":
            return b"changed blob" if problem == "git_blob" else source.read_bytes()
        return "wrong" if problem == "revision" else module.COMMIT
    monkeypatch.setattr(module.subprocess, "check_output", check_output)
    def git_diff(*args, **kwargs):
        if problem == "dirty":
            raise subprocess.CalledProcessError(1, args[0])
    monkeypatch.setattr(module.subprocess, "run", git_diff)
    if problem == "file":
        Path(records[0]["path"]).write_text("changed")
    if problem == "submission":
        extra.write_text("changed")
    return report


@pytest.mark.parametrize("problem", [None, "admission_pin", "readback_pin", "source_pin", "revision", "dirty",
    "file", "submission", "readback_semantics", "git_revision", "git_blob", "admission_binding", "source_binding", "conflict"])
def test_source_gate(tmp_path, monkeypatch, problem):
    report = source_fixture(tmp_path, monkeypatch, problem)
    if problem:
        with pytest.raises((ValueError, subprocess.CalledProcessError)):
            module.verify_candidates(tmp_path)
    else:
        verified = module.verify_candidates(tmp_path)
        assert verified["admission"] == report and verified["arm"]["label"] == "cpm_high"


def allocation(monkeypatch):
    for key, value in dict(SLURM_JOB_ID="99999", SLURM_CPUS_PER_TASK="32",
        SLURM_MEM_PER_NODE="196608", SLURM_JOB_NODELIST="bizon").items():
        monkeypatch.setenv(key, value)
    monkeypatch.delenv("SLURM_ARRAY_TASK_ID", raising=False)
    monkeypatch.setattr(module, "PYTHON", sys.executable)


@pytest.mark.parametrize("key,value", [("SLURM_JOB_ID", ""), ("SLURM_CPUS_PER_TASK", "2"),
    ("SLURM_MEM_PER_NODE", "65536"), ("SLURM_JOB_NODELIST", "dgx"), ("SLURM_ARRAY_TASK_ID", "1")])
def test_driver_allocation_precedes_inputs(tmp_path, monkeypatch, key, value):
    allocation(monkeypatch)
    monkeypatch.setenv(key, value)
    with pytest.raises(ValueError, match="standalone"):
        module.run(tmp_path, "unused")


def test_helper_interpreter_not_authorized(tmp_path, monkeypatch):
    allocation(monkeypatch)
    monkeypatch.setattr(module, "PYTHON", str(tmp_path / "different-python"))
    with pytest.raises(ValueError, match="full runtime"):
        module.run(tmp_path, "unused")


def test_baseline_runtime_failure_is_not_bypassed(tmp_path, monkeypatch):
    from benchmark_tools import run_qfo_parameter_phylogeny as baseline_module

    verified = {"admission": "candidate admission passed"}
    monkeypatch.setattr(module, "verify_candidates", lambda *args: verified)
    def reject(*args):
        raise ValueError("Package inventory changed: orthohmm")
    monkeypatch.setattr(baseline_module, "verify_baseline", reject)
    with pytest.raises(ValueError, match="Package inventory changed"):
        module.verify_sources(tmp_path)
    assert not (tmp_path / module.OUTPUT).exists()


@pytest.mark.parametrize("outcome", ["success", "fresh", "admission_failure", "failed", "exception", "changed", "protocol"])
def test_execution_contract(tmp_path, monkeypatch, outcome):
    from benchmark_tools import run_qfo_factorial_cell as command_module
    from benchmark_tools import run_simulation_methods as execution_module

    allocation(monkeypatch)
    launcher = tmp_path / "launcher"
    launcher.mkdir()
    protocol = tmp_path / module.PROTOCOL
    protocol.parent.mkdir(parents=True)
    protocol.write_text("prospective protocol")
    original = next(cell for cell in plan_cells(tmp_path / "baseline", tmp_path / "fastas",
        tmp_path / "executor/benchmark_tools/replay_phylogeny.py", 32) if cell["label"] == "p1_c1_r1")
    report = report_fixture(tmp_path)
    report.update(preparation={"sha256": "preparation-pin"}, protocol={"sha256": "admission-protocol-pin"})
    arm = module.select_arm(tmp_path, report)
    before = copy.deepcopy(original)
    verified = {"launcher": str(launcher), "prepared": str(tmp_path / "executor"), "original": original,
        "arm": arm, "admission": report, "admission_executor": str(tmp_path / "validator"),
        "environment": {}, "checked_records": [],
        "manifest": {"input_fastas": [], "environment_overrides": {"OMP_NUM_THREADS": "1"}}}
    checks, launches, admissions = [], [], []
    def verify(*args):
        checks.append(args)
        return {} if outcome == "changed" and len(checks) > 1 else verified
    monkeypatch.setattr(module, "verify_sources", verify)
    monkeypatch.setattr(command_module, "native_command", lambda cell, *args: (cell["argv"], []))
    monkeypatch.setattr(execution_module, "execution_environment", lambda *args: ({}, {}))
    def admission_run(command, **kwargs):
        admissions.append(command)
        assert command[2] == str(tmp_path / "validator/benchmark_tools/admit_helper_cpm_candidates.py")
        assert command[command.index("--preparation-sha256") + 1] == "preparation-pin"
        assert command[command.index("--protocol-sha256") + 1] == "admission-protocol-pin"
        if outcome == "admission_failure":
            raise subprocess.CalledProcessError(1, command)
        Path(command[-1]).write_text(json.dumps({} if outcome == "fresh" else report))
    monkeypatch.setattr(module.subprocess, "run", admission_run)
    def execute(dataset, order, env, output, inputs, provenance):
        launches.append(dataset)
        assert Path.cwd() == launcher and order == ["candidate_cpm_high"]
        argv = dataset["methods"][order[0]]["argv"]
        assert argv[argv.index("--candidate-clusters") + 1] == arm["partition"]["path"]
        assert argv[argv.index("--membership-constraints") + 1] == arm["constraints"]["path"]
        assert argv[argv.index("--species-tree-mode") + 1] == "infer"
        assert argv[argv.index("--checkpoint-source") + 1] == original["argv"][original["argv"].index("--output-directory") + 1]
        assert "--species-tree" not in argv and "--official-benchmark" not in argv
        assert env == {"OMP_NUM_THREADS": "1", "PYTHONPATH": str(launcher)}
        assert "incremental shared-host" in provenance["scope"]
        if outcome == "exception":
            raise RuntimeError("native interruption")
        return {"failed_methods": order if outcome == "failed" else []}
    monkeypatch.setattr(execution_module, "execute", execute)
    cwd = Path.cwd()
    sha = "wrong" if outcome == "protocol" else record(protocol)["sha256"]
    if outcome == "success":
        result = module.run(tmp_path, sha)
        assert result["status"] == "complete_pending_native_validation"
    else:
        with pytest.raises((ValueError, RuntimeError, subprocess.CalledProcessError)):
            module.run(tmp_path, sha)
    destination = tmp_path / module.OUTPUT
    if outcome == "protocol":
        assert not destination.exists() and not checks and not admissions and not launches
        return
    result = json.loads((destination / "postflight.json").read_text())
    assert result["status"] == ("complete_pending_native_validation" if outcome == "success" else "failed")
    assert all(result[key] is False for key in ("accuracy_evaluated", "native_outputs_validated", "publication_ready"))
    assert len(admissions) == 1 and len(launches) == (0 if outcome in ("fresh", "admission_failure") else 1)
    assert Path.cwd() == cwd and original == before
    with pytest.raises(FileExistsError):
        module.run(tmp_path, sha)


def test_no_scientific_import_before_gate():
    code = ("import sys; import benchmark_tools.run_helper_cpm_phylogeny; "
            "assert not any(n == 'orthohmm' or n.startswith('orthohmm.') for n in sys.modules)")
    subprocess.run([sys.executable, "-S", "-B", "-c", code], check=True)


def test_launch_policy():
    text = (Path(module.__file__).parent / "results/qfo_cpm_helper_phylogeny_20261001.sh").read_text()
    for fragment in ("--nodelist=bizon", "--cpus-per-task=32", "--mem=192G", "--time=24:00:00",
        "--no-requeue", "unset PYTHONNOUSERSITE", module.PYTHON, "--protocol-sha256", "diff --quiet HEAD"):
        assert fragment in text
    assert "--array" not in text and "--dependency" not in text and "dgx" not in text
