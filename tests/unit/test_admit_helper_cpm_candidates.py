import copy
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import admit_helper_cpm_candidates as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


ACCOUNTING = "JobID|JobIDRaw|State|ExitCode|Elapsed|AllocCPUS|ReqMem|NodeList\n22385|22385|COMPLETED|0:0|00:10:00|2|64G|bizon\n"


@pytest.mark.parametrize("problem", [None, "missing", "duplicate", "running", "failed", "exit", "raw", "cpus", "memory", "host"])
def test_completion_contract(problem):
    text = ACCOUNTING
    if problem == "missing":
        text = text.splitlines()[0] + "\n"
    elif problem == "duplicate":
        text += text.splitlines()[1] + "\n"
    elif problem in ("running", "failed"):
        text = text.replace("COMPLETED", problem.upper())
    elif problem == "exit":
        text = text.replace("|0:0|", "|1:0|")
    elif problem == "raw":
        text = text.replace("22385|22385", "22385|other")
    elif problem == "cpus":
        text = text.replace("|2|64G|", "|32|64G|")
    elif problem == "memory":
        text = text.replace("64G", "32G")
    elif problem == "host":
        text = text.replace("bizon", "other")
    if problem:
        with pytest.raises(ValueError, match="completed"):
            module.completed(text)
    else:
        assert module.completed(text)["JobID"] == "22385"


@pytest.mark.parametrize("problem", [None, "status", "arm", "job", "executor", "source", "builder", "context",
    "inputs", "handoff", "accuracy", "publication", "label", "delta", "parameters", "calls", "report",
    "expansion", "nan", "negative", "bool_time"])
def test_preparation_contract(problem):
    parameters = dict(min_norm=.03, min_margin=1.5)
    source, builder, inputs, context = {"source": True}, {"builder": True}, [{"input": True}], {"context": True}
    expansion = {"parameters": parameters}
    report = dict(status="recovered_cpm_candidates_prepared_pending_admission", arm="cpm_high", job_id=module.JOB,
        executor_commit=module.COMMIT, source=source, builder=builder, inputs=inputs, context=context,
        seed_handoff="explicit_helper_runtime_seed_amendment", accuracy_evaluated=False, publication_ready=False,
        candidate_parameter_control=dict(label="control", delta={}, applied_parameters=parameters,
            engine_calls=1, engine_fixed_profile_report=expansion), candidate_arm={"expansion": expansion}, incremental_seconds=1.)
    report = copy.deepcopy(report)
    control = report["candidate_parameter_control"]
    mutations = dict(status=(report, "status", "running"), arm=(report, "arm", "cpm_low"),
        job=(report, "job_id", "other"), executor=(report, "executor_commit", "other"), source=(report, "source", {}),
        builder=(report, "builder", {}), context=(report, "context", {}), inputs=(report, "inputs", []),
        handoff=(report, "seed_handoff", "historical_admission_22155"), accuracy=(report, "accuracy_evaluated", True),
        publication=(report, "publication_ready", True), label=(control, "label", "norm_low"),
        delta=(control, "delta", {"min_norm": .024}), parameters=(control, "applied_parameters", {}),
        calls=(control, "engine_calls", True), report=(control, "engine_fixed_profile_report", {"parameters": {}}),
        expansion=(report["candidate_arm"], "expansion", {}), nan=(report, "incremental_seconds", float("nan")),
        negative=(report, "incremental_seconds", -1.), bool_time=(report, "incremental_seconds", True))
    if problem:
        target, field, value = mutations[problem]
        target[field] = value
        with pytest.raises(ValueError):
            module.validate_report(report, {"JobIDRaw": module.JOB}, context, parameters, source, builder, inputs)
    else:
        module.validate_report(report, {"JobIDRaw": module.JOB}, context, parameters, source, builder, inputs)


@pytest.mark.parametrize("state", ["PENDING", "RUNNING", "FAILED"])
def test_no_file_or_scientific_access_before_completion(tmp_path, monkeypatch, state):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: ACCOUNTING.replace("COMPLETED", state))
    def forbidden(*args):
        raise AssertionError("No scientific import before actual completion")
    monkeypatch.setattr(module, "scientific_imports", forbidden)
    output = tmp_path / "admission.json"
    with pytest.raises(ValueError, match="completed"):
        module.admit(tmp_path, "fixture", "fixture", output)
    assert not output.exists()


@pytest.mark.parametrize("problem", [None, "submission_sha", "status", "job", "executor", "revision", "flag",
                                    "source", "gate", "changed_record"])
def test_submission_pins(tmp_path, monkeypatch, problem):
    executor = tmp_path / module.EXECUTOR
    def file(path, value="fixture"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(value)
        return record(path)
    source = file(executor / "benchmark_tools/prepare_recovered_cpm_candidates.py")
    gate = file(executor / "benchmark_tools/helper_recovered_cpm_seed_evidence.py")
    readback = file(tmp_path / "readback.json")
    failure = file(tmp_path / "failure.json")
    data = dict(status="source_corrected_explicit_helper_seed_candidate_job_submitted", job_id=module.JOB,
        executor=str(executor), executor_commit=module.COMMIT, seed_handoff="explicit_helper_runtime_seed_amendment",
        candidate_admitted=False, accuracy_evaluated=False, publication_ready=False, previous_failed_attempt=failure,
        seed_readback=readback, source_records=[source, gate])
    if problem in ("status", "job", "executor", "flag"):
        key, value = dict(status=("status", "submitted"), job=("job_id", "other"),
                          executor=("executor", "other"), flag=("candidate_admitted", True))[problem]
        data[key] = value
    item = file(tmp_path / module.SUBMISSION, json.dumps(data))
    monkeypatch.setattr(module, "SUBMISSION_SHA", "wrong" if problem == "submission_sha" else item["sha256"])
    monkeypatch.setattr(module, "SOURCE_SHA", "wrong" if problem == "source" else source["sha256"])
    monkeypatch.setattr(module, "GATE_SHA", "wrong" if problem == "gate" else gate["sha256"])
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "other" if problem == "revision" else module.COMMIT)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    if problem == "changed_record":
        Path(readback["path"]).write_text("changed")
    if problem:
        with pytest.raises(ValueError):
            module.pinned_submission(tmp_path)
    else:
        actual, frozen, observed_source, observed_gate, observed_executor, checked = module.pinned_submission(tmp_path)
        assert actual == data and frozen == item and observed_executor == executor
        assert observed_source == source and observed_gate == gate and failure in checked


@pytest.mark.parametrize("problem", [None, "wrong_origin"])
def test_scientific_selection_keeps_origin_guard(tmp_path, monkeypatch, problem):
    paths = [tmp_path / "orthohmm/orthohmm.py", tmp_path / "orthohmm/accuracy.py"]
    if problem:
        paths[0] = tmp_path / "wrong/orthohmm.py"
    imports = iter(SimpleNamespace(__file__=str(path)) for path in paths)
    monkeypatch.setattr(module.importlib, "import_module", lambda name: next(imports))
    before = list(sys.path)
    if problem:
        with pytest.raises(ValueError, match="scientific import"):
            module.scientific_imports(tmp_path)
    else:
        assert Path(module.scientific_imports(tmp_path).__file__) == paths[-1]
    assert sys.path == before


@pytest.mark.parametrize("different", [False, True])
def test_exact_seed_report_reconstruction_uses_frozen_module_path(tmp_path, monkeypatch, different):
    path = tmp_path / "benchmark_tools/helper_recovered_cpm_seed_evidence.py"
    path.parent.mkdir()
    path.write_text("fixture")
    selected = dict(protocol={"sha256": "protocol"}, checked_records=[{"path": str(path)}])
    calls = []
    def fresh(root, readback, digest, protocol):
        calls.append((root, readback, digest, protocol))
        return {} if different else selected
    loaded = SimpleNamespace(READBACK="readback.json", READBACK_SHA="fixed", evidence=fresh)
    specs = []
    spec = SimpleNamespace(loader=SimpleNamespace(exec_module=lambda actual: specs.append(actual)))
    monkeypatch.setattr(module.importlib.util, "spec_from_file_location",
                        lambda name, origin: spec if origin == str(path) else None)
    monkeypatch.setattr(module.importlib.util, "module_from_spec", lambda actual: loaded)
    if different:
        with pytest.raises(ValueError, match="handoff differs"):
            module.replay_seed_evidence(tmp_path, tmp_path, record(path), {"recovery_admission": selected})
    else:
        assert module.replay_seed_evidence(tmp_path, tmp_path, record(path), {"recovery_admission": selected}) == selected
    assert specs == [loaded]
    assert calls == [(tmp_path, tmp_path / "readback.json", "fixed", "protocol")]


@pytest.mark.parametrize("problem", [None, "preparation_sha", "protocol_sha", "report", "runtime", "numeric",
                                    "auditor", "universe", "content", "merge", "changed_input", "seed_after", "runtime_after"])
def test_independent_admission_orchestration(tmp_path, monkeypatch, problem):
    from benchmark_tools import audit_accuracy_checkpoint as numeric_module
    from benchmark_tools import audit_candidate_arm as content_module
    from benchmark_tools import checked_replay_payload_worker as corrected
    from benchmark_tools import cpm_replay_context as context_module
    from benchmark_tools import run_simulation_methods as frozen
    from benchmark_tools import trace_ob_families as trace
    from benchmark_tools import verify_qfo_replay_launcher as runtime_module
    monkeypatch.setattr(module, "GENES", 4)
    def file(relative, value="fixture"):
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(value)
        return record(path)
    protocol = file(module.PROTOCOL)
    executor = tmp_path / module.EXECUTOR
    source = file(module.EXECUTOR + "/benchmark_tools/prepare_recovered_cpm_candidates.py")
    gate = file(module.EXECUTOR + "/benchmark_tools/helper_recovered_cpm_seed_evidence.py")
    builder = file(module.EXECUTOR + "/benchmark_tools/prepare_qfo_cpm_candidates.py")
    auditor = file(module.EXECUTOR + "/benchmark_tools/audit_accuracy_checkpoint.py")
    decoder = file("decoder.py")
    names = file("names.txt", "a\nb\nc\nd\n" if problem != "universe" else "a\nb\nb\nd\n")
    seed = file("seed.txt", "a b\nc d\n")
    checkpoint = file("checkpoint/manifest.json")
    plan_record = file("plan.json")
    baseline_record = file("benchmarks/results/qfo_corrected_factorial_v1/manifest.json")
    runtime_record = file("benchmark_tools/results/publication_native_runtime_20260916.json")
    submission_record = file(module.SUBMISSION)
    file(f"benchmarks/work/qfo_cpm_helper_candidates_{module.JOB}.log")
    file(f"benchmarks/work/qfo_cpm_helper_candidates_{module.JOB}.time.txt")
    parameters = dict(min_norm=.03, min_margin=1.5)
    context = dict(checked_records=[])
    plan = dict(runtime={"fixed": True}, checkpoint_manifest=checkpoint)
    numeric = dict(status="numeric_checkpoint_verified", manifest=checkpoint, summary={"genes": 4, "species": 78},
                   accuracy_evaluated=False, auditor=auditor)
    baseline = dict(runtime_before=plan["runtime"], runtime_after=plan["runtime"], numeric_checkpoint=numeric,
        candidate_arms={"p1_c1": {"expansion": {"parameters": parameters}}})
    admitted = dict(seed_partition=seed, checked_records=[seed], report_record=file("seed_admission.json"),
                    readback=file("seed_readback.json"))
    partition = file(str(Path(module.PREPARATION).parent / "candidate/partition.txt"), "a b c d\n")
    events = file(str(Path(module.PREPARATION).parent / "candidate/merges.json"), "[]")
    expansion = dict(parameters=parameters)
    arm = dict(candidate_partition=partition, membership_constraints=events,
               output_files=[partition, events], expansion=expansion)
    content = dict(status="fixture", checked_records=[seed, partition, events])
    helpers = [record(path) for path in sorted((executor / "benchmark_tools").glob("*.py"))]
    inputs = [plan_record, baseline_record, runtime_record, builder, seed, *helpers]
    report = dict(status="recovered_cpm_candidates_prepared_pending_admission", arm="cpm_high", job_id=module.JOB,
        executor_commit=module.COMMIT, source=source, builder=builder, inputs=inputs, context=context,
        seed_handoff="explicit_helper_runtime_seed_amendment", accuracy_evaluated=False, publication_ready=False,
        candidate_parameter_control=dict(label="control", delta={}, applied_parameters=parameters,
            engine_calls=1, engine_fixed_profile_report=expansion), candidate_arm=arm, incremental_seconds=1.,
        runtime_before=plan["runtime"], runtime_after=plan["runtime"], numeric_checkpoint=copy.deepcopy(numeric),
        content_audit=content)
    if problem == "report":
        report["seed_handoff"] = "wrong"
    elif problem == "auditor":
        report["numeric_checkpoint"]["auditor"] = {}
    elif problem == "content":
        report["content_audit"] = {}
    preparation = file(module.PREPARATION, json.dumps(report))
    monkeypatch.setattr(module, "pinned_submission", lambda *a: ({}, submission_record, source, gate, executor, [submission_record]))
    seed_calls = []
    def fresh(*args):
        seed_calls.append(args)
        return {} if problem == "seed_after" and len(seed_calls) == 2 else admitted
    monkeypatch.setattr(module, "replay_seed_evidence", fresh)
    monkeypatch.setattr(module, "scientific_imports", lambda *a: SimpleNamespace(__file__=decoder["path"]))
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: ACCOUNTING)
    monkeypatch.setattr(corrected, "corrected_evidence", lambda *a: (plan, plan_record, {}, names))
    monkeypatch.setattr(context_module, "evidence", lambda *a: context)
    monkeypatch.setattr(frozen, "read_frozen", lambda *a: baseline)
    monkeypatch.setattr(numeric_module, "audit", lambda *a: {} if problem == "numeric" else numeric)
    runtime_calls = []
    def runtime(*args):
        runtime_calls.append(args)
        return {} if problem == "runtime" or problem == "runtime_after" and len(runtime_calls) == 2 else plan["runtime"]
    monkeypatch.setattr(runtime_module, "verify", runtime)
    monkeypatch.setattr(content_module, "audit", lambda *a: content)
    monkeypatch.setattr(trace, "partition", lambda path, *a: ({"0": set("ab"), "1": set("cd")} if path == Path(seed["path"])
                                                             else {"0": set("abcd")}, {}))
    reconstructions = []
    def reconstruction(*args):
        reconstructions.append(args)
        if problem == "merge":
            raise ValueError("Invalid merge reconstruction")
        if problem == "changed_input":
            Path(seed["path"]).write_text("changed")
    monkeypatch.setattr(trace, "validate_merge_reconstruction", reconstruction)
    output = tmp_path / "admission.json"
    arguments = (tmp_path, "wrong" if problem == "preparation_sha" else preparation["sha256"],
                 "wrong" if problem == "protocol_sha" else protocol["sha256"], output)
    if problem:
        with pytest.raises((ValueError, KeyError)):
            module.admit(*arguments)
        assert not output.exists()
    else:
        result = module.admit(*arguments)
        assert result == json.loads(output.read_text())
        assert result["status"] == "cpm_helper_recovered_candidates_admitted_unscored"
        assert result["candidate_admitted"] is True
        assert result["accuracy_evaluated"] is result["downstream_admitted"] is result["publication_ready"] is False
        assert result["verification"] == dict(genes=4, seed_groups=2, candidate_groups=1, reconstructed_merges=0)
        assert len(reconstructions) == 1 and len(seed_calls) == len(runtime_calls) == 2
        with pytest.raises(FileExistsError):
            module.admit(*arguments)


def test_cli_help_imports_no_scientific_package(tmp_path):
    root = Path(module.__file__).resolve().parent.parent
    script = "import runpy,sys; sys.argv=[sys.argv[1],\"--help\"];\ntry: runpy.run_path(sys.argv[0],run_name=\"__main__\")\nexcept SystemExit as error: assert error.code==0\nassert not any(name==\"orthohmm\" or name.startswith(\"orthohmm.\") for name in sys.modules)"
    done = subprocess.run([sys.executable, "-I", "-B", "-c", script,
        str(root / "benchmark_tools/admit_helper_cpm_candidates.py")], cwd=tmp_path, capture_output=True, text=True)
    assert done.returncode == 0, done.stdout + done.stderr


def test_launch_script_matches_prospective_admission_protocol():
    root = Path(module.__file__).resolve().parent.parent
    path = root / "benchmark_tools/results/qfo_cpm_helper_candidate_admission_20260930.sh"
    done = subprocess.run(["bash", "-n", str(path)], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    content = path.read_text()
    assert "--protocol-sha256 " + record(root / module.PROTOCOL)["sha256"] in content
    assert "--preparation-sha256 2e62321bf7a6d410f8ddafff783aaea979f1fec753fd6415b2b439ee7e2f54e6" in content
    for flag in ("--nodelist=bizon", "--ntasks=1", "--cpus-per-task=2", "--mem=64G", "--time=01:00:00", "--no-requeue"):
        assert "#SBATCH " + flag in content
    assert "/home/bizon/anaconda3/bin/python -B" in content
    assert "--dependency" not in content
    assert "dgx" not in content.lower()
