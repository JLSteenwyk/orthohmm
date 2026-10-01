import copy
import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import run_private_helper_cpm_phylogeny as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.prepare_orthobench_factorial import plan_cells
from tests.unit.test_readback_qfo_private_phylogeny_control import accounting
from tests.unit.test_run_helper_cpm_phylogeny import report_fixture


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return record(path)


def authority_fixture(root, monkeypatch, problem):
    from benchmark_tools.readback_qfo_private_phylogeny_control import completed

    executor = root / module.EXECUTOR
    source = executor / "benchmark_tools/admit_qfo_private_phylogeny_control.py"
    source.parent.mkdir(parents=True)
    source.write_text("frozen admission source")
    probe = root / "benchmark_tools/readback_qfo_private_phylogeny_control.py"
    probe.parent.mkdir(parents=True)
    probe.write_text("frozen readback source")
    data = root / "bound-evidence"
    data.write_text("unchanged evidence")
    submission = dict(job_id="22389", executor=str(executor), executor_commit=module.ADMISSION_COMMIT,
        source_records=[record(source)], parent_submission=record(data), preserved_failed_attempt=record(data))
    submission_pin = write_json(root / "submission.json", submission)
    report = dict(status="private_qfo_phylogeny_deployment_admitted_unscored", job_id="22387",
        scheduler=completed(accounting())[0], native_group_integrity=dict(native_outputs_validated=True),
        recovered_cpm_inference_authorized=True, accuracy_evaluated=False, scoring_admitted=False,
        controlled_timing=False, publication_ready=False, native_pair_count=5959560,
        native_comparison={}, cache_use={}, source=record(source), checked_records=[record(data)])
    admission_pin = write_json(root / module.ADMISSION, report)
    readback = dict(status="private_qfo_phylogeny_admission_independently_read_back_unscored",
        admission=admission_pin, source=record(probe), source_revision=module.ADMISSION_COMMIT,
        scheduler=completed(accounting()), native_pair_count=5959560, native_comparison={}, cache_use={},
        recovered_cpm_inference_authorized=True, accuracy_evaluated=False, scoring_admitted=False,
        controlled_timing=False, publication_ready=False, submission=submission_pin,
        source_git_bindings=[{**record(source), "git_revision": module.ADMISSION_COMMIT}])
    changes = {"readback_status": (readback, "status", "pending"), "readback_source": (readback, "source", record(data)),
        "readback_revision": (readback, "source_revision", "wrong"), "readback_scheduler": (readback, "scheduler", []),
        "readback_pairs": (readback, "native_pair_count", 0), "readback_comparison": (readback, "native_comparison", {"x": 0}),
        "readback_cache": (readback, "cache_use", {"x": 0}), "readback_authority": (readback, "recovered_cpm_inference_authorized", 1),
        "readback_scoring": (readback, "scoring_admitted", True), "report_status": (report, "status", "pending"),
        "report_authority": (report, "recovered_cpm_inference_authorized", False), "report_accuracy": (report, "accuracy_evaluated", True),
        "report_groups": (report, "native_group_integrity", dict(native_outputs_validated=False)),
        "report_source": (report, "source", record(data)), "git_binding": (readback, "source_git_bindings", [])}
    if problem in changes:
        target, key, value = changes[problem]
        target[key] = value
    admission_pin = write_json(root / module.ADMISSION, report)
    readback["admission"] = admission_pin
    readback_pin = write_json(root / module.READBACK, readback)
    monkeypatch.setattr(module, "ADMISSION_SHA", "wrong" if problem == "admission_sha" else admission_pin["sha256"])
    monkeypatch.setattr(module, "READBACK_SHA", "wrong" if problem == "readback_sha" else readback_pin["sha256"])
    monkeypatch.setattr(module, "READBACK_SOURCE_SHA", "wrong" if problem == "probe_sha" else record(probe)["sha256"])
    def command(argv, **kwargs):
        if argv[0] == "sacct":
            return accounting().replace("COMPLETED", "RUNNING") if problem == "active" else accounting()
        if argv[-2] == "show":
            return b"wrong blob" if problem == "git_blob" else source.read_bytes()
        return "wrong" if problem == "revision" else module.ADMISSION_COMMIT
    monkeypatch.setattr(module.subprocess, "check_output", command)
    def diff(*args, **kwargs):
        if problem == "dirty":
            raise subprocess.CalledProcessError(1, args[0])
    monkeypatch.setattr(module.subprocess, "run", diff)
    if problem == "file":
        data.write_text("changed")


@pytest.mark.parametrize("problem", [None, "probe_sha", "active", "admission_sha", "readback_sha", "readback_status",
    "readback_source", "readback_revision", "readback_scheduler", "readback_pairs", "readback_comparison", "readback_cache",
    "readback_authority", "readback_scoring", "report_status", "report_authority", "report_accuracy", "report_groups",
    "report_source", "git_binding", "git_blob", "revision", "dirty", "file"])
def test_private_control_authority(tmp_path, monkeypatch, problem):
    authority_fixture(tmp_path, monkeypatch, problem)
    if problem:
        with pytest.raises((ValueError, subprocess.CalledProcessError)):
            module.verify_private_control(tmp_path)
    else:
        assert module.verify_private_control(tmp_path)["scheduler"][1]["JobID"] == "22389"


def test_original_candidate_and_private_baseline_gates_reused(tmp_path, monkeypatch):
    from benchmark_tools import run_helper_cpm_phylogeny as candidates
    from benchmark_tools import qfo_private_phylogeny_environment as baseline

    monkeypatch.setattr(module, "verify_private_control", lambda *args: {"checked_records": ["private"]})
    monkeypatch.setattr(candidates, "verify_candidates", lambda *args: {"checked_records": ["candidate"]})
    monkeypatch.setattr(baseline, "verify_baseline", lambda *args: {"checked_records": ["baseline"]})
    assert module.verify_sources(tmp_path)["checked_records"] == ["private", "candidate", "baseline"]
    def reject(*args):
        raise ValueError("private runtime changed")
    monkeypatch.setattr(baseline, "verify_baseline", reject)
    with pytest.raises(ValueError, match="runtime changed"):
        module.verify_sources(tmp_path)


def context(root):
    original = next(cell for cell in plan_cells(root / "baseline", root / "fastas",
        root / "executor/benchmark_tools/replay_phylogeny.py", 32) if cell["label"] == "p1_c1_r1")
    report = report_fixture(root)
    report.update(preparation={"sha256": "preparation-pin"}, protocol={"sha256": "admission-protocol-pin"})
    from benchmark_tools.run_helper_cpm_phylogeny import select_arm
    return dict(candidates=dict(arm=select_arm(root, report), admission=report, admission_executor=str(root / "validator")),
        baseline=dict(original=original, launcher=str(root / "launcher"), prepared=str(root / "executor"),
            environment={"tool_entrypoints": {"orthohmm_python": {"absolute_path": str(root / "private-python")}}},
            manifest=dict(input_fastas=[], environment_overrides={"OMP_NUM_THREADS": "1"})), checked_records=[])


@pytest.mark.parametrize("problem", [None, "cpu", "mode", "supplied", "scoring", "checkpoint"])
def test_private_command_only_changes_interpreter(tmp_path, monkeypatch, problem):
    from benchmark_tools import run_qfo_factorial_cell as frozen
    verified = context(tmp_path)
    original = copy.deepcopy(verified)
    def command(cell, *args):
        argv = list(cell["argv"])
        if problem in ("cpu", "mode", "checkpoint"):
            flag = {"cpu": "--cpu", "mode": "--species-tree-mode", "checkpoint": "--checkpoint-source"}[problem]
            argv[argv.index(flag) + 1] = "wrong"
        elif problem:
            argv.extend(["--species-tree" if problem == "supplied" else "--official-benchmark", "wrong"])
        return argv, []
    monkeypatch.setattr(frozen, "native_command", command)
    if problem:
        with pytest.raises(ValueError):
            module.native_command(verified, tmp_path / "new")
    else:
        cell, planned, argv, _ = module.native_command(verified, tmp_path / "new")
        assert argv[0] == str(tmp_path / "private-python") and argv[1:] == planned["argv"][1:]
        assert cell["argv"] == argv and verified == original


def allocation(monkeypatch):
    for key, value in dict(SLURM_JOB_ID="99999", SLURM_CPUS_PER_TASK="32", SLURM_MEM_PER_NODE="196608",
                          SLURM_JOB_NODELIST="bizon").items():
        monkeypatch.setenv(key, value)
    monkeypatch.delenv("SLURM_ARRAY_TASK_ID", raising=False)
    monkeypatch.setattr(module, "PYTHON", sys.executable)


@pytest.mark.parametrize("key,value", [("SLURM_JOB_ID", ""), ("SLURM_CPUS_PER_TASK", "2"),
    ("SLURM_MEM_PER_NODE", "65536"), ("SLURM_JOB_NODELIST", "dgx"), ("SLURM_ARRAY_TASK_ID", "1")])
def test_bad_allocation_precedes_inputs(tmp_path, monkeypatch, key, value):
    allocation(monkeypatch)
    monkeypatch.setenv(key, value)
    with pytest.raises(ValueError, match="standalone"):
        module.run(tmp_path, "unused")


@pytest.mark.parametrize("outcome", ["success", "fresh", "admission_failed", "lookup_failed", "failed", "exception", "changed", "protocol"])
def test_execution_contract(tmp_path, monkeypatch, outcome):
    from benchmark_tools import run_simulation_methods as execution
    from benchmark_tools import qfo_private_phylogeny_environment as private
    from benchmark_tools import run_qfo_factorial_cell as frozen
    allocation(monkeypatch)
    verified = context(tmp_path)
    Path(verified["baseline"]["launcher"]).mkdir()
    protocol = tmp_path / module.PROTOCOL
    protocol.parent.mkdir(parents=True)
    protocol.write_text("prospective private protocol")
    checks, launches = [], []
    def verify(*args):
        checks.append(args)
        return {} if outcome == "changed" and len(checks) > 1 else verified
    monkeypatch.setattr(module, "verify_sources", verify)
    monkeypatch.setattr(frozen, "native_command", lambda cell, *args: (list(cell["argv"]), []))
    monkeypatch.setattr(execution, "execution_environment", lambda *args: ({"LD_PRELOAD": "unsafe", "PYTHONHOME": "wrong"}, {}))
    def lookup(*args):
        if outcome == "lookup_failed":
            raise ValueError("lookup failed")
        return {"checked_records": []}
    monkeypatch.setattr(private, "inspect_launcher", lookup)
    def candidate_admission(argv, **kwargs):
        assert argv[0] == sys.executable and "admit_helper_cpm_candidates.py" in argv[2]
        if outcome == "admission_failed":
            raise subprocess.CalledProcessError(1, argv)
        Path(argv[-1]).write_text(json.dumps({} if outcome == "fresh" else verified["candidates"]["admission"]))
    monkeypatch.setattr(module.subprocess, "run", candidate_admission)
    def execute(dataset, order, env, *args):
        launches.append(dataset)
        assert order == ["candidate_cpm_high"] and Path.cwd() == Path(verified["baseline"]["launcher"])
        assert dataset["methods"][order[0]]["argv"][0] == str(tmp_path / "private-python")
        assert env["PYTHONNOUSERSITE"] == env["PYTHONDONTWRITEBYTECODE"] == "1"
        assert env["NUMBA_CACHE_DIR"] == str(tmp_path / module.OUTPUT / "numba_cache")
        assert "LD_PRELOAD" not in env and "PYTHONHOME" not in env
        if outcome == "exception":
            raise RuntimeError("interrupted")
        return {"failed_methods": order if outcome == "failed" else []}
    monkeypatch.setattr(execution, "execute", execute)
    cwd = Path.cwd()
    sha = "wrong" if outcome == "protocol" else record(protocol)["sha256"]
    if outcome == "success":
        assert module.run(tmp_path, sha)["status"] == "complete_pending_native_validation"
    else:
        with pytest.raises((ValueError, RuntimeError, subprocess.CalledProcessError)):
            module.run(tmp_path, sha)
    output = tmp_path / module.OUTPUT
    if outcome == "protocol":
        assert not output.exists() and not checks and not launches
    else:
        post = json.loads((output / "postflight.json").read_bytes())
        assert post["status"] == ("complete_pending_native_validation" if outcome == "success" else "failed")
        assert all(post[k] is False for k in ("accuracy_evaluated", "native_outputs_validated", "publication_ready"))
        assert len(launches) == (0 if outcome in ("fresh", "admission_failed", "lookup_failed") else 1)
        with pytest.raises(FileExistsError):
            module.run(tmp_path, sha)
    assert Path.cwd() == cwd


def test_cli_imports_no_scientific_code():
    code = ("import sys; import benchmark_tools.run_private_helper_cpm_phylogeny; "
            "assert not any(n == 'orthohmm' or n.startswith('orthohmm.') for n in sys.modules)")
    subprocess.run([sys.executable, "-S", "-B", "-c", code], check=True)


def test_frozen_standalone_script():
    path = Path(module.__file__).parent / "results/qfo_private_cpm_phylogeny_20261001.sh"
    subprocess.run(["bash", "-n", str(path)], check=True)
    text = path.read_text()
    for token in ("--cpus-per-task=32", "--mem=192G", "--time=24:00:00", "--no-requeue",
                  "--nodelist=bizon", "--protocol-sha256", "run_private_helper_cpm_phylogeny.py"):
        assert token in text
    assert "--dependency" not in text and "--array" not in text and "dgx" not in text
