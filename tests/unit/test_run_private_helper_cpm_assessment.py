import copy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import run_private_helper_cpm_assessment as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit import test_prepare_private_helper_cpm_pairs as conversion_tests


JOB = "24002"


def fixture(tmp_path, monkeypatch):
    data = conversion_tests.setup(tmp_path, monkeypatch)
    root = Path(module.__file__).resolve().parents[1]
    executor = tmp_path / module.CONVERTER_EXECUTOR
    names = ("benchmark_tools/prepare_private_helper_cpm_pairs.py", module.converter.PROTOCOL,
        "benchmark_tools/results/qfo_private_cpm_pair_conversion_20261001.sh",
        "tests/unit/test_prepare_private_helper_cpm_pairs.py")
    for name in names:
        path = executor / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes((root / name).read_bytes())
    data.protocol.write_bytes((root / module.converter.PROTOCOL).read_bytes())
    monkeypatch.setattr(module.converter, "__file__", str(executor / names[0]))
    environment_path = tmp_path / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    infrastructure = tmp_path / "runtime.cfg"
    infrastructure.write_text("# pinned fixture\n")
    environment = dict(status="local_qfo_assessment_environment_frozen", accuracy_evaluated=False,
        reference_files=[record(data.mapping)], source=record(infrastructure),
        execution_config=record(infrastructure), singularity_config=record(infrastructure),
        pipeline_files=[], java_files=[], images=[], executables=[], singularity_support=[],
        pipeline=str(tmp_path / "qfo_benchmark/benchmark-webservice"),
        environment_overrides={"JAVA_HOME": "/frozen/java", "NXF_OFFLINE": "true"})
    environment_path.write_text(json.dumps(environment))
    monkeypatch.setattr(module.converter, "ENV_SHA", record(environment_path)["sha256"])
    stage = conversion_tests.invoke(tmp_path, data)
    stage_path = data.output / "results.json"
    submitted = dict(status="private_recovered_qfo_high_cpm_pair_conversion_submitted", job_id=JOB,
        executor=str(executor), executor_commit=module.CONVERTER_COMMIT, executor_clean=True,
        native_inference=False, accuracy_evaluated=False, controlled_timing=False, publication_ready=False,
        native_admission=stage["native_admission"], parent_submission=record(data.submission),
        source_records=[record(executor / name) for name in names],
        submission_argv=["sbatch", "--parsable", str(executor / names[2]), str(executor), module.CONVERTER_COMMIT,
            module.CONVERTER_PROTOCOL_SHA, stage["admission_job"], stage["native_admission"]["sha256"],
            record(data.submission)["sha256"]])
    submission = tmp_path / f"benchmark_tools/results/qfo_private_cpm_pair_conversion_submission_{JOB}.json"
    submission.write_text(json.dumps(submitted))
    protocol = tmp_path / module.PROTOCOL
    protocol.write_text("# Prospective assessment fixture\n")
    monkeypatch.setattr(module.converter, "accounting", lambda job: conversion_tests.accounting().replace(conversion_tests.JOB, job))
    monkeypatch.setattr(module.subprocess, "check_output", lambda argv, **kw:
        module.CONVERTER_COMMIT if module.CONVERTER_EXECUTOR in argv[2] else module.converter.ADMITTER_COMMIT)
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(returncode=0))
    monkeypatch.setattr(module, "command_for", lambda root, stage, env, work, results:
        ["frozen-nextflow", "--challenges_ids", "GO EC VGNC SwissTrees TreeFam-A FAS", "--input", stage["filtered_pairs"]["path"]])
    return SimpleNamespace(data=data, stage=stage, stage_path=stage_path, submitted=submitted,
        submission=submission, protocol=protocol, environment=environment, environment_path=environment_path,
        output=tmp_path / module.OUTPUT, work=tmp_path / module.WORK, results=tmp_path / module.RESULTS)


def verify(root, data):
    return module.verify_conversion(root, JOB, record(data.stage_path)["sha256"], record(data.submission)["sha256"])


def prepare(root, data):
    return module.prepare(root, JOB, record(data.stage_path)["sha256"],
                          record(data.submission)["sha256"], record(data.protocol)["sha256"])


def test_full_verified_conversion_and_fresh_assessment_preparation(tmp_path, monkeypatch):
    data = fixture(tmp_path, monkeypatch)
    value = verify(tmp_path, data)
    assert value["stage"] == data.stage
    assert value["scheduler"]["JobIDRaw"] == JOB
    assert data.stage["pairs"] in value["checked_records"]
    prepared, again = prepare(tmp_path, data)
    assert again == value
    assert prepared["arm"] == "cpm_high" and prepared["index"] == 1
    assert prepared["conversion_job"] == JOB and prepared["converter_commit"] == module.CONVERTER_COMMIT
    assert prepared["stage"]["semantics"] == "native phylogenetically inferred pairs"
    assert prepared["accuracy_admitted"] is False and prepared["publication_ready"] is False
    assert prepared["environment_manifest"] == record(data.environment_path)
    assert prepared["command"][2] == "GO EC VGNC SwissTrees TreeFam-A FAS"
    assert all(not path.exists() for path in (data.output, data.work, data.results))


@pytest.mark.parametrize("key,value", [("status", "old_cpm_status"), ("arm", "cpm_low"),
    ("index", True), ("participant", "other"), ("semantics", "RootHOG cliques"), ("job_id", "999"),
    ("source", {}), ("accuracy_evaluated", True), ("scoring_admitted", 0), ("publication_ready", 0),
    ("controlled_timing", True), ("written_pairs", 3), ("retained_pairs", True),
    ("removed_mapping_pairs", 1), ("filtered_pairs", {}), ("checked_records", [])])
def test_incomplete_or_promoted_stage_rejected(tmp_path, monkeypatch, key, value):
    data = fixture(tmp_path, monkeypatch)
    stage = copy.deepcopy(data.stage)
    stage[key] = value
    with pytest.raises((ValueError, KeyError)):
        module.validate_stage(stage, {"JobIDRaw": JOB}, data.stage["source"])


def test_active_conversion_precedes_any_file_reads(tmp_path, monkeypatch):
    monkeypatch.setattr(module.converter, "accounting", lambda job:
        conversion_tests.accounting(state="RUNNING").replace(conversion_tests.JOB, job))
    monkeypatch.setattr(module, "read_frozen", lambda *args: pytest.fail("Read live conversion"))
    with pytest.raises(ValueError, match="completed"):
        module.verify_conversion(tmp_path, JOB, "", "")


@pytest.mark.parametrize("key,value", [("status", "wrong"), ("job_id", "999"),
    ("executor_commit", "wrong"), ("executor_clean", False), ("source_records", []),
    ("submission_argv", []), ("native_admission", {}), ("parent_submission", {}),
    ("native_inference", True), ("accuracy_evaluated", 0)])
def test_conversion_submission_rejection(tmp_path, monkeypatch, key, value):
    data = fixture(tmp_path, monkeypatch)
    data.submitted[key] = value
    data.submission.write_text(json.dumps(data.submitted))
    with pytest.raises(ValueError):
        verify(tmp_path, data)


@pytest.mark.parametrize("problem", ["preflight", "fresh", "count_sidecar", "native_command",
    "native_context", "path", "missing_native_submission", "mapping", "protocol", "conversion_pin", "dirty"])
def test_changed_context_counts_preflight_or_mapping_rejected(tmp_path, monkeypatch, problem):
    data = fixture(tmp_path, monkeypatch)
    if problem == "preflight":
        path = data.data.output / "preflight.json"
        obj = json.loads(path.read_text())
        obj["status"] = "running"
        path.write_text(json.dumps(obj))
        old = module.bound_record(data.stage["checked_records"], path)
        data.stage["checked_records"][data.stage["checked_records"].index(old)] = record(path)
    elif problem == "fresh":
        path = data.data.output / "native_admission_recheck.json"
        path.write_text("{}")
        old = data.stage["native_admission_recheck"]
        data.stage["native_admission_recheck"] = record(path)
        data.stage["checked_records"][data.stage["checked_records"].index(old)] = record(path)
    elif problem == "count_sidecar":
        path = Path(data.stage["conversion_counts"]["path"])
        path.write_text("{}")
        old = data.stage["conversion_counts"]
        data.stage["conversion_counts"] = record(path)
        data.stage["checked_records"][data.stage["checked_records"].index(old)] = record(path)
    elif problem == "native_command":
        data.stage["native_recheck_command"] = ["wrong"]
    elif problem == "native_context":
        data.stage["candidate_admission"] = {}
    elif problem == "path":
        data.stage["pairs"] = {**data.stage["pairs"], "path": "/alternate/pairs.tsv"}
    elif problem == "missing_native_submission":
        ref = module.bound_record(data.stage["checked_records"], data.data.submission)
        data.stage["checked_records"].remove(ref)
    elif problem == "mapping":
        data.environment["reference_files"] = []
        data.environment_path.write_text(json.dumps(data.environment))
        monkeypatch.setattr(module.converter, "ENV_SHA", record(data.environment_path)["sha256"])
    elif problem == "protocol":
        data.protocol.write_text("changed")
    elif problem == "conversion_pin":
        with pytest.raises(ValueError):
            module.verify_conversion(tmp_path, JOB, "0" * 64, record(data.submission)["sha256"])
        return
    else:
        def run(*args, **kwargs):
            raise module.subprocess.CalledProcessError(1, args[0])
        monkeypatch.setattr(module.subprocess, "run", run)
    data.stage_path.write_text(json.dumps(data.stage))
    with pytest.raises((ValueError, module.subprocess.CalledProcessError)):
        if problem == "protocol":
            module.prepare(tmp_path, JOB, record(data.stage_path)["sha256"], record(data.submission)["sha256"], "0" * 64)
        else:
            prepare(tmp_path, data)
    assert not data.output.exists()


def test_missing_or_conflicting_binding_rejected():
    ref = dict(path="/record", bytes=1, sha256="a")
    assert module.bound_record([ref, ref], Path("/record")) == ref
    for records in ([], [ref, {**ref, "sha256": "b"}]):
        with pytest.raises(ValueError):
            module.bound_record(records, Path("/record"))


@pytest.mark.parametrize("relative", [module.OUTPUT, module.WORK, module.RESULTS])
@pytest.mark.parametrize("symlink", [False, True])
def test_namespace_gate_before_expensive_checks(tmp_path, monkeypatch, relative, symlink):
    path = tmp_path / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    if symlink:
        path.symlink_to(tmp_path / "absent")
    else:
        path.mkdir()
    monkeypatch.setattr(module, "verify_conversion", lambda *a: pytest.fail("Repeated full verification after collision"))
    with pytest.raises(FileExistsError):
        module.prepare(tmp_path, JOB, "", "", "")


@pytest.mark.parametrize("outcome", ["success", "exit", "exception", "source_changed", "context_changed", "inventory_error"])
def test_execution_retains_failure_without_admitting_scores(tmp_path, monkeypatch, outcome):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(module.sys, "executable", module.converter.PYTHON)
    for key, value in {"SLURM_JOB_ID": "24003", "SLURM_CPUS_PER_TASK": "8",
                      "SLURM_JOB_NODELIST": "bizon", "SLURM_MEM_PER_NODE": "65536"}.items():
        monkeypatch.setenv(key, value)
    monkeypatch.delenv("SLURM_ARRAY_TASK_ID", raising=False)
    source = tmp_path / "source.py"
    source.write_text("# unchanged\n")
    output, results = tmp_path / "execution", tmp_path / "scores"
    results.mkdir()
    (results / "score.json").write_text("{}")
    prepared = dict(status="prepared_unrun", cwd=str(output), results=str(results), source=record(source),
        verified_records=[], command=["frozen-nextflow", "six-endpoints"], environment_overrides={"NXF_OFFLINE": "true"},
        accuracy_admitted=False, controlled_timing=False, publication_ready=False)
    verified = {"fixed": "context"}
    monkeypatch.setattr(module, "prepare", lambda *a: (prepared.copy(), verified))
    monkeypatch.setattr(module, "verify_conversion", lambda *a: {} if outcome == "context_changed" else verified)
    def execute(command, **kwargs):
        assert command == prepared["command"] and kwargs["cwd"] == output
        assert kwargs["env"]["NXF_OFFLINE"] == "true"
        assert kwargs["env"]["NXF_SINGULARITY_CACHEDIR"].endswith("qfo_benchmark/scoring/container_cache")
        if outcome == "exception":
            raise RuntimeError("Interrupted")
        if outcome == "source_changed":
            source.write_text("# changed\n")
        return SimpleNamespace(returncode=1 if outcome in ("exit", "inventory_error") else 0)
    monkeypatch.setattr(module.subprocess, "run", execute)
    if outcome == "inventory_error":
        original = module.record
        def record_error(path):
            if Path(path).name == "score.json":
                raise OSError("Unavailable partial artifact")
            return original(path)
        monkeypatch.setattr(module, "record", record_error)
    if outcome == "success":
        module.run(tmp_path, JOB, "converted", "submitted", "protocol")
    else:
        with pytest.raises((RuntimeError, ValueError)):
            module.run(tmp_path, JOB, "converted", "submitted", "protocol")
    saved = json.loads((output / "results.json").read_text())
    assert saved["status"] == ("process_succeeded_pending_independent_admission" if outcome == "success" else "failed")
    assert saved["accuracy_admitted"] is False and saved["publication_ready"] is False
    assert (output / "preflight.json").exists() and (output / "scoring.log").exists()
    if outcome == "inventory_error":
        assert saved["error_type"] == "RuntimeError" and saved["output_capture_error"]["type"] == "OSError"


@pytest.mark.parametrize("key,value", [("SLURM_JOB_ID", ""), ("SLURM_CPUS_PER_TASK", "2"),
    ("SLURM_JOB_NODELIST", "other"), ("SLURM_MEM_PER_NODE", "131072"), ("SLURM_ARRAY_TASK_ID", "1")])
def test_allocation_before_preparation(tmp_path, monkeypatch, key, value):
    for name, setting in {"SLURM_JOB_ID": "24003", "SLURM_CPUS_PER_TASK": "8",
                         "SLURM_JOB_NODELIST": "bizon", "SLURM_MEM_PER_NODE": "65536"}.items():
        monkeypatch.setenv(name, setting)
    monkeypatch.setenv(key, value)
    monkeypatch.setattr(module.sys, "executable", module.converter.PYTHON)
    monkeypatch.setattr(module, "prepare", lambda *a: pytest.fail("Prepared without allocation"))
    with pytest.raises(ValueError, match="controller"):
        module.run(tmp_path, JOB, "", "", "")


def test_original_command_builder_preserves_endpoints_without_resume():
    from benchmark_tools.run_qfo_recovered_assessment import command_for

    env = dict(execution_config=dict(path="/config"), pipeline="/pipeline")
    stage = dict(filtered_pairs=dict(path="/pairs.qfo.tsv"), participant="ohmm_qfo_parameter_cpm_high")
    argv = command_for(Path("/project"), stage, env, Path("/w/qcpp1"), Path("/results"))
    assert argv[argv.index("--challenges_ids") + 1] == "GO EC VGNC SwissTrees TreeFam-A FAS"
    assert argv[argv.index("--event_year") + 1] == "2020"
    assert argv[argv.index("--input") + 1] == "/pairs.qfo.tsv"
    assert "-resume" not in argv and "--challenges_ids" in argv


def test_mapping_mismatch_after_valid_conversion_rejected(tmp_path, monkeypatch):
    data = fixture(tmp_path, monkeypatch)
    verified = verify(tmp_path, data)
    monkeypatch.setattr(module, "verify_conversion", lambda *args: verified)
    data.environment["reference_files"] = []
    data.environment_path.write_text(json.dumps(data.environment))
    monkeypatch.setattr(module.converter, "ENV_SHA", record(data.environment_path)["sha256"])
    with pytest.raises(ValueError, match="mapping differ"):
        prepare(tmp_path, data)
    assert not data.output.exists()


def test_wrong_converter_revision_or_source_rejected(tmp_path, monkeypatch):
    data = fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(module.subprocess, "check_output", lambda *args, **kw: "wrong")
    with pytest.raises(ValueError, match="revision"):
        verify(tmp_path, data)
    monkeypatch.setattr(module, "CONVERTER_SHA", "0" * 64)
    with pytest.raises(ValueError, match="source changed"):
        verify(tmp_path, data)


def test_wrong_assessment_controller_or_cwd_rejected(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    for name, setting in {"SLURM_JOB_ID": "24003", "SLURM_CPUS_PER_TASK": "8",
                         "SLURM_JOB_NODELIST": "bizon", "SLURM_MEM_PER_NODE": "65536"}.items():
        monkeypatch.setenv(name, setting)
    monkeypatch.delenv("SLURM_ARRAY_TASK_ID", raising=False)
    monkeypatch.setattr(module.sys, "executable", "/unreviewed/python")
    with pytest.raises(ValueError, match="controller"):
        module.run(tmp_path, JOB, "", "", "")
    monkeypatch.setattr(module.sys, "executable", module.converter.PYTHON)
    with pytest.raises(ValueError, match="verification directory"):
        module.run(tmp_path / "other", JOB, "", "", "")


def test_script_requires_exact_resources_and_input_bindings():
    script = Path(module.__file__).parent / "results/qfo_private_cpm_assessment_20261001.sh"
    text = script.read_text()
    for required in ("--nodelist=bizon", "--cpus-per-task=8", "--mem=64G", "--time=24:00:00",
        "--no-requeue", "run_private_helper_cpm_assessment.py", "--conversion-job", "--conversion-sha256",
        "--conversion-submission-sha256", "--protocol-sha256", "/home/bizon/anaconda3/bin/python -B"):
        assert required in text
    for forbidden in ("--array", "--dependency", "--gres", "--exclusive", "sbatch ", "while ", "until "):
        assert forbidden not in text
