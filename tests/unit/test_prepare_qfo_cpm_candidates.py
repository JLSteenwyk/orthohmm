import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import prepare_qfo_cpm_candidates as module
from benchmark_tools.audit_candidate_arm import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_prepare_qfo_candidate_neighborhood import fixture


def test_variant_seed_and_unchanged_candidate_parameters(tmp_path):
    engine, original, seen, baseline, output, _ = fixture(tmp_path)
    variant = tmp_path / "variant.txt"
    variant.write_text("a\nb c\n")
    result = module.build(engine, baseline["expansion"]["parameters"], record(variant),
                          list("abc"), [0, 1, 2], (), output / "candidate", audit)
    assert result["candidate_arm"]["seed_partition"] == record(variant)
    assert Path(result["candidate_arm"]["candidate_partition"]["path"]).read_bytes() == variant.read_bytes()
    assert result["content_audit"]["genes"] == 3
    assert result["candidate_parameter_control"]["delta"] == {}
    assert result["candidate_parameter_control"]["applied_parameters"] == baseline["expansion"]["parameters"]
    assert len(seen) == 1 and engine.merge_supported_satellite_candidate_clusters is original
    with pytest.raises(FileExistsError):
        module.build(engine, baseline["expansion"]["parameters"], record(variant),
                     list("abc"), [0, 1, 2], (), output / "candidate", audit)


def test_engine_failure_restores_function_and_preserves_seed(tmp_path):
    engine, _, _, baseline, output, _ = fixture(tmp_path)
    def fail(*args, **kwargs):
        raise RuntimeError("native failure")
    engine.merge_supported_satellite_candidate_clusters = fail
    with pytest.raises(RuntimeError, match="native failure"):
        module.build(engine, baseline["expansion"]["parameters"], baseline["seed_partition"],
                     list("abc"), [0, 1, 2], (), output / "candidate", audit)
    assert engine.merge_supported_satellite_candidate_clusters is fail
    assert (output / "candidate/orthohmm_working_res/orthohmm_edges_clustered.txt").read_text() == "a b\nc\n"


@pytest.mark.parametrize("state", ["RUNNING", "PENDING", "FAILED"])
def test_scheduler_before_missing_outputs(tmp_path, monkeypatch, state):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        f"JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n21962_0|21963|{state}|0:0|00:01:00|bizon|2\n")
    with pytest.raises(ValueError, match="completed"):
        module.replay_evidence(tmp_path, 0, {})


@pytest.mark.parametrize("index", [-1, 2, True, "0"])
def test_bad_index(tmp_path, index):
    with pytest.raises(ValueError, match="index"):
        module.replay_evidence(tmp_path, index, {})


@pytest.mark.parametrize("problem", [None, "source", "revision", "arm", "context", "seed", "stages", "changed_seed", "cpus"])
def test_replay_evidence(tmp_path, monkeypatch, problem):
    executor = tmp_path / "benchmarks/work/publication_qfo_cpm_variant_admission_v1"
    source = executor / "benchmark_tools/admit_qfo_cpm_variant.py"
    source.parent.mkdir(parents=True)
    source.write_text("fixture")
    monkeypatch.setattr(module, "ADMISSION_SHA", "wrong" if problem == "source" else record(source)["sha256"])
    output = tmp_path / "variant"
    seed = output / "replay/orthogroups_profiles_refined.txt"
    seed.parent.mkdir(parents=True)
    seed.write_text("a b\nc\n")
    context = {"output_root": str(output), "expected_stages": ["multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"]}
    report = {"status": "cpm_variant_replay_admitted_unscored", "index": 0,
        "arm": "control" if problem == "arm" else "cpm_low", "source": record(source),
        "context": {} if problem == "context" else context,
        "accuracy_evaluated": False, "publication_ready": False, "checked_records": [record(seed)],
        "coverage": [{"label": label, "output": record(seed)} for label in context["expected_stages"]]}
    if problem == "seed":
        report["coverage"][-1]["output"]["path"] = str(tmp_path / "baseline.txt")
    elif problem == "stages":
        report["coverage"].reverse()
    path = tmp_path / "benchmarks/work/qfo_cpm_variant_admission_21962_0.json"
    path.write_text(json.dumps(report))
    if problem == "changed_seed":
        seed.write_text("changed")
    def check_output(command, **kwargs):
        if command[0] == "sacct":
            cpus = "32" if problem == "cpus" else "2"
            return f"JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n21962_0|21963|COMPLETED|0:0|00:01:00|bizon|{cpus}\n"
        return "changed" if problem == "revision" else module.ADMISSION_COMMIT
    monkeypatch.setattr(module.subprocess, "check_output", check_output)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    if problem:
        with pytest.raises(ValueError):
            module.replay_evidence(tmp_path, 0, context)
    else:
        result = module.replay_evidence(tmp_path, 0, context)
        assert result["seed_partition"] == record(seed)
        assert result["report"] == report


def test_unscheduled_preparation_rejected(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="scheduled"):
        module.prepare(tmp_path, 0)


def test_numeric_auditor_cannot_preload_executor_scientific_package(tmp_path):
    code = r'''
import os
from pathlib import Path
import sys
sys.path.insert(0, sys.argv[1])
from benchmark_tools import prepare_qfo_cpm_candidates as m
from benchmark_tools import checked_replay_payload_worker as payload
from benchmark_tools import cpm_replay_context as context
from benchmark_tools import verify_qfo_replay_launcher as runtime
root = Path(sys.argv[2])
os.environ['SLURM_ARRAY_TASK_ID'] = '0'
m.require_environment = lambda env: None
payload.corrected_evidence = lambda *a: ({'runtime': {}}, {}, {}, {})
context.evidence = lambda *a: {'cwd': str(root / 'launcher')}
m.replay_evidence = lambda *a: {}
m.read_frozen = lambda *a: {'candidate_arms': {'p1_c1': {'expansion': {'parameters': {}}}},
                          'runtime_before': {}, 'runtime_after': {}}
runtime.verify = lambda *a: {}
original = m.importlib.import_module
def load(name, *args, **kwargs):
    if name == 'orthohmm.orthohmm':
        assert 'orthohmm' not in sys.modules, 'scientific package preloaded before launcher selection'
        assert sys.path[0] == str(root / 'launcher')
        raise RuntimeError('frozen-import-boundary-reached')
    return original(name, *args, **kwargs)
m.importlib.import_module = load
try:
    m.prepare(root, 0)
except RuntimeError as error:
    assert str(error) == 'frozen-import-boundary-reached'
else:
    raise AssertionError('Scientific import boundary not exercised')
'''
    subprocess.run([sys.executable, "-I", "-c", code, str(Path(module.__file__).resolve().parent.parent),
                    str(tmp_path)], check=True)


@pytest.mark.parametrize("problem", [None, "fresh", "engine", "runtime", "checkpoint", "import"])
def test_preparation_orchestration(tmp_path, monkeypatch, problem):
    def write(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))
        return record(path)
    monkeypatch.setattr(module, "require_environment", lambda *a: None)
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "0")
    monkeypatch.setenv("SLURM_JOB_ID", "fixture")
    launcher = tmp_path / "launcher"
    checkpoint = write(tmp_path / "checkpoint/manifest.json", {})
    numeric = {"summary": {"genes": 984137, "species": 78}}
    params = {"min_norm": .03, "min_margin": 1.5}
    baseline = {"runtime_before": {}, "runtime_after": {}, "numeric_checkpoint": numeric,
                "candidate_arms": {"p1_c1": {"expansion": {"parameters": params}}}}
    baseline_record = write(tmp_path / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json", baseline)
    monkeypatch.setattr(module, "BASELINE_SHA", baseline_record["sha256"])
    plan = {"runtime": {}, "checkpoint_manifest": checkpoint}
    plan_record = write(tmp_path / "plan.json", plan)
    runtime_path = tmp_path / "benchmark_tools/results/publication_native_runtime_20260916.json"
    write(runtime_path, {})
    context = {"cwd": str(launcher), "checked_records": []}
    seed = write(tmp_path / "variant_seed.txt", {})
    admitted = {"executor": str(tmp_path / "validator"), "report": {"fixture": True},
                "seed_partition": seed, "checked_records": [seed]}
    monkeypatch.setattr("benchmark_tools.checked_replay_payload_worker.corrected_evidence", lambda *a: (plan, plan_record, {}, {}))
    monkeypatch.setattr("benchmark_tools.cpm_replay_context.evidence", lambda *a: context)
    monkeypatch.setattr(module, "replay_evidence", lambda *a: admitted)
    runtime_calls, audit_calls = [], []
    def verify(*args):
        runtime_calls.append(args)
        return {"changed": True} if problem == "runtime" and len(runtime_calls) > 1 else {}
    def audit(*args):
        audit_calls.append(args)
        return {} if problem == "checkpoint" and len(audit_calls) > 1 else numeric
    monkeypatch.setattr("benchmark_tools.verify_qfo_replay_launcher.verify", verify)
    monkeypatch.setattr("benchmark_tools.audit_accuracy_checkpoint.audit", audit)
    names = [f"g{i}" for i in range(984137)]
    engine = SimpleNamespace(__file__=str(launcher / "orthohmm/orthohmm.py"))
    accuracy = SimpleNamespace(__file__=str(launcher / "orthohmm/accuracy.py"),
                               load_accuracy_checkpoint=lambda *a, **k: (names, [], [], [], []))
    if problem == "import":
        engine.__file__ = str(tmp_path / "wrong/orthohmm.py")
    original_import = module.importlib.import_module
    def load(name, *args, **kwargs):
        if name == "orthohmm.orthohmm":
            return engine
        if name == "orthohmm.accuracy":
            return accuracy
        return original_import(name, *args, **kwargs)
    monkeypatch.setattr(module.importlib, "import_module", load)
    builds = []
    def build(*args):
        builds.append(args)
        assert args[1] == params and args[2] == seed and len(args[3]) == 984137
        if problem == "engine":
            raise RuntimeError("candidate failure")
        item = write(args[6] / "partition.txt", {})
        return {"candidate_arm": {"output_files": [item]}, "content_audit": {"fixture": True}}
    monkeypatch.setattr(module, "build", build)
    def run(command, **kwargs):
        write(Path(command[-1]), {} if problem == "fresh" else admitted["report"])
    monkeypatch.setattr(module.subprocess, "run", run)
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "frozen")
    destination = tmp_path / "benchmarks/results/qfo_cpm_candidates_v1/cpm_low"
    if problem:
        with pytest.raises((ValueError, RuntimeError)):
            module.prepare(tmp_path, 0)
        if problem == "import":
            assert not destination.exists()
        else:
            report = json.loads((destination / "manifest.json").read_text())
            assert report["status"] == "failed"
            if problem == "fresh":
                assert not builds
    else:
        report = module.prepare(tmp_path, 0)
        assert report["status"] == "cpm_candidates_prepared_pending_admission"
        assert report["accuracy_evaluated"] is False and report["publication_ready"] is False
        assert len(builds) == 1
        assert json.loads((destination / "manifest.json").read_text()) == report
        with pytest.raises(FileExistsError):
            module.prepare(tmp_path, 0)
