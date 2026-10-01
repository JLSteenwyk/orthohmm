import copy
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import admit_qfo_private_phylogeny_control as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def accounting(state="COMPLETED", exit_code="0:0", cpus="32", memory="192G", node="bizon", step="COMPLETED"):
    return ("JobID|JobIDRaw|State|ExitCode|Elapsed|AllocCPUS|ReqMem|NodeList\n"
        f"22387|22387|{state}|{exit_code}|00:20:00|{cpus}|{memory}|{node}\n"
        f"22387.batch|22387.batch|{step}|0:0|00:20:00|32|192G|bizon\n")


@pytest.mark.parametrize("values", [{}, {"state": "RUNNING"}, {"state": "FAILED"}, {"state": "CANCELLED"},
    {"exit_code": "1:0"}, {"cpus": "2"}, {"memory": "64G"}, {"node": "dgx"}, {"step": "RUNNING"}, {"step": "FAILED"}])
def test_scheduler_gate(values):
    if values:
        with pytest.raises(ValueError):
            module.completed(accounting(**values))
    else:
        assert module.completed(accounting())["JobIDRaw"] == "22387"


@pytest.mark.parametrize("text", ["", accounting().split("22387.batch")[0],
    accounting() + accounting().splitlines()[1] + "\n"])
def test_missing_or_ambiguous_scheduler(text):
    with pytest.raises(ValueError):
        module.completed(text)


def test_active_parent_precedes_data_reads(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "accounting", lambda: accounting(state="RUNNING"))
    with pytest.raises(ValueError, match="completed"):
        module.admit(tmp_path, "unused", "unused", tmp_path / "admission.json")
    assert not (tmp_path / "admission.json").exists()


def identity_fixture():
    expected = {"job_id": "22387"}
    inputs = {"status": "ready", "inputs": []}
    post = {"status": "private_qfo_baseline_parity_complete_pending_admission", "native_comparison": {},
        "accuracy_evaluated": False, "publication_ready": False,
        "recovered_cpm_inference_authorized": False, "controlled_timing": False}
    status = {"provenance": expected, "verified_inputs": inputs, "dataset": module.LABEL,
        "methods": {module.LABEL: {}}, "failed_methods": [], "status": "finished_pending_native_validation",
        "accuracy_evaluated": False, "native_outputs_validated": False}
    return expected, inputs, post, status


@pytest.mark.parametrize("problem", [None, "preflight", "provenance", "inputs", "dataset", "methods", "failed", "status", "accuracy", "native", "post_status", "post_flag", "post_error"])
def test_complete_execution_identity(problem):
    expected, inputs, post, status = identity_fixture()
    preflight = copy.deepcopy(expected)
    mutations = {"preflight": (preflight, "job_id", "wrong"), "provenance": (status, "provenance", {}),
        "inputs": (status, "verified_inputs", {}), "dataset": (status, "dataset", "other"),
        "methods": (status, "methods", {}), "failed": (status, "failed_methods", [module.LABEL]),
        "status": (status, "status", "running"), "accuracy": (status, "accuracy_evaluated", 0),
        "native": (status, "native_outputs_validated", True), "post_status": (post, "status", "failed"),
        "post_flag": (post, "recovered_cpm_inference_authorized", True), "post_error": (post, "error", "failure")}
    if problem:
        target, key, value = mutations[problem]
        target[key] = value
        with pytest.raises(ValueError):
            module.execution_identity(preflight, post, status, expected, inputs)
    else:
        module.execution_identity(preflight, post, status, expected, inputs)


@pytest.mark.parametrize("problem", [None, "bool", "negative", "candidate", "total", "hits", "remapped", "species"])
def test_checkpoint_accounting(problem):
    summary = dict(candidate_families=351739, reconciled_families=100, bypassed_families=351639,
        checkpoint_hits=90, remapped_checkpoint_hits=5, species_tree_families=10, species_tree_checkpoint_hit=True)
    changes = {"bool": ("checkpoint_hits", True), "negative": ("remapped_checkpoint_hits", -1),
        "candidate": ("candidate_families", 351738), "total": ("bypassed_families", 0),
        "hits": ("checkpoint_hits", 101), "remapped": ("remapped_checkpoint_hits", 91),
        "species": ("species_tree_checkpoint_hit", 1)}
    if problem:
        key, value = changes[problem]
        summary[key] = value
        with pytest.raises(ValueError):
            module.cache_summary(summary)
    else:
        assert module.cache_summary(summary) == summary


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return record(path)


def admission_fixture(root, monkeypatch, problem):
    from benchmark_tools import run_simulation_methods as execution
    from benchmark_tools import validate_simulation_outputs as process_module
    from benchmark_tools import validate_factorial_native as native_module
    from benchmark_tools import admit_qfo_corrected_factorial_cell as ownership_module
    from benchmark_tools import admit_qfo_factorial_cell as pairs_module
    from benchmark_tools.inspect_native_python_lookup import PROBE, scientific_origins
    from benchmark_tools.run_qfo_private_phylogeny_control import compare_outputs, NATIVE_FILES

    executor = root / module.EXECUTOR
    helpers = []
    for name in ("run_qfo_private_phylogeny_control.py", "qfo_private_phylogeny_environment.py"):
        source = executor / "benchmark_tools" / name
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_bytes((Path(module.__file__).parent / name).read_bytes())
        helpers.append(record(source))
    helpers.sort(key=lambda p: p["path"])
    control_relative = "benchmark_tools/results/QFO_PRIVATE_PHYLOGENY_CONTROL_PROTOCOL_20261001.md"
    control_protocol = root / control_relative
    control_protocol.parent.mkdir(parents=True)
    control_protocol.write_text("control protocol")
    executor_protocol = executor / control_relative
    executor_protocol.parent.mkdir(parents=True)
    executor_protocol.write_bytes(control_protocol.read_bytes())
    monkeypatch.setattr(module, "CONTROL_PROTOCOL_SHA", record(control_protocol)["sha256"])
    protocol = root / module.PROTOCOL
    protocol.parent.mkdir(parents=True, exist_ok=True)
    protocol.write_text("admission protocol")
    source = record(executor / "benchmark_tools/run_qfo_private_phylogeny_control.py")
    submission = {"job_id": module.JOB, "status": "private_full_qfo_phylogeny_control_submitted",
        "executor_commit": module.COMMIT, "executor": str(executor), "source_records": [*helpers, record(executor_protocol)],
        "accuracy_evaluated": False, "recovered_cpm_inference_authorized": False, "controlled_timing": False,
        "publication_ready": False, "deployment": record(protocol), "original_native_admission": record(protocol),
        "preserved_shared_rejection": record(protocol)}
    submission_pin = write_json(root / module.SUBMISSION, submission)
    output = root / module.OUTPUT
    launcher = root / "launcher"
    prepared = root / "prepared"
    python = root / "private-python"
    python.write_text("private executable")
    for path in (launcher / "orthohmm/phylogeny_pipeline.py", launcher / "benchmark_tools/replay_phylogeny.py"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("frozen source")
    equivalence = []
    for name in ("replay_phylogeny.py", "orthobench_diagnostics.py"):
        executed = launcher / "benchmark_tools" / name
        executed.write_text("frozen source")
        original_source = prepared / "benchmark_tools" / name
        original_source.parent.mkdir(parents=True, exist_ok=True)
        original_source.write_bytes(executed.read_bytes())
        equivalence.append({"prepared": record(original_source), "executed": record(executed)})
    old = root / "baseline/orthohmm_phylogeny"
    new = output / "output/orthohmm_phylogeny"
    old.mkdir(parents=True)
    new.mkdir(parents=True)
    for name in NATIVE_FILES:
        (old / name).write_text(name)
        (new / name).write_text(name + ("changed" if problem == "parity" and name == NATIVE_FILES[0] else ""))
    candidate = root / "candidate.txt"
    candidate.write_text("a b\n")
    argv = [str(python), str(launcher / "benchmark_tools/replay_phylogeny.py"), "--output-directory", str(output / "output"),
        "--json", str(output / "metrics.json"), "--checkpoint-source", str(old.parent)]
    original = {"label": "p1_c1_r1", "candidate_partition": str(candidate)}
    verified = {"original": original, "launcher": str(launcher), "prepared": str(prepared),
        "environment": {}, "manifest": {"input_fastas": []}, "checked_records": []}
    calls = []
    def verify(*args):
        calls.append(args)
        return {} if problem == "runtime_changed" and len(calls) > 1 else verified
    env_module = SimpleNamespace(verify_baseline=verify, control_command=lambda *args: (argv, equivalence))
    monkeypatch.setattr(module, "control_environment", lambda *args: (env_module, verify()))
    monkeypatch.setattr(module, "accounting", lambda: accounting())
    monkeypatch.setattr(module.subprocess, "check_output", lambda *args, **kwargs: "wrong" if problem == "revision" else module.COMMIT)
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: None)
    monkeypatch.setattr(execution, "execution_environment", lambda *args: ({}, {}))
    preflight = {"source": source, "helpers": helpers, "protocol": record(control_protocol), "verified": verified,
        "executed_argv": argv, "launcher_equivalence": equivalence, "resolved_tools": {}, "cwd": str(launcher), "job_id": module.JOB,
        "scope": "Entire frozen QfO p1_c1_r1 private deployment parity; validated checkpoint reuse; unscored incremental shared-host control"}
    _, inputs, post, status = identity_fixture()
    if problem == "protocol_path":
        preflight["protocol"] = record(executor_protocol)
    status.update(provenance=preflight, verified_inputs=inputs)
    comparisons = compare_outputs(old, new)
    post["native_comparison"] = comparisons
    if problem == "comparison_report":
        post["native_comparison"] = {}
    if problem == "execution":
        status["status"] = "running"
    write_json(output / "preflight.json", preflight)
    write_json(output / "postflight.json", post)
    write_json(output / "execution/status.json", status)
    report = {"modules": {"orthohmm.phylogeny_pipeline": str(launcher / "orthohmm/phylogeny_pipeline.py"),
        "benchmark_tools.replay_phylogeny": str(launcher / "benchmark_tools/replay_phylogeny.py")},
        "mapped_files": [], "executable": str(python), "cwd": str(launcher), "dont_write_bytecode": True,
        "pycache_prefix": str(output / "bytecode_cache"),
        "requested": ["orthohmm.phylogeny_pipeline", "benchmark_tools.replay_phylogeny"]}
    lookup = {"report": report, "origin": scientific_origins(report, "orthohmm", launcher / "orthohmm"),
        "checked_records": [record(path) for path in sorted({*report["modules"].values(), str(python)})],
        "continuous_enforcement": False}
    lookup_command = [str(python), "-B", "-c", PROBE, json.dumps(report["requested"])]
    if problem == "lookup":
        lookup["continuous_enforcement"] = True
    if problem == "lookup_request":
        report["requested"] = []
    write_json(output / "lookup.json", lookup)
    write_json(output / "lookup_process.json", dict(returncode=0, command=lookup_command, stdout=json.dumps(report), stderr=""))
    summary = dict(candidate_families=351739, reconciled_families=100, bypassed_families=351639,
        checkpoint_hits=90, remapped_checkpoint_hits=5, species_tree_families=10, species_tree_checkpoint_hit=True, ortholog_pairs=1)
    write_json(new / "provenance_manifest.json", {"results": summary})
    write_json(output / "metrics.json", {})
    monkeypatch.setattr(process_module, "verify_process", lambda *args: {str(new / NATIVE_FILES[0])})
    monkeypatch.setattr(native_module, "validate_native_cell", lambda *args, **kwargs: dict(
        native_manifest=record(new / "provenance_manifest.json"), native_metrics=record(output / "metrics.json"),
        species_tree=record(new / "species_tree.rooted.nwk"), native_outputs_validated=True))
    monkeypatch.setattr(ownership_module, "gene_ownership", lambda *args: ({"a": "s1", "b": "s2"}, {"a": 0, "b": 0}))
    monkeypatch.setattr(pairs_module, "check_pairs", lambda *args: 0 if problem == "no_pairs" else 1)
    return submission_pin["sha256"], record(protocol)["sha256"], calls


@pytest.mark.parametrize("problem", [None, "revision", "execution", "protocol_path", "lookup", "lookup_request", "parity", "comparison_report", "no_pairs", "runtime_changed"])
def test_end_to_end_readonly_gate(tmp_path, monkeypatch, problem):
    submission_sha, protocol_sha, calls = admission_fixture(tmp_path, monkeypatch, problem)
    destination = tmp_path / "admission.json"
    if problem:
        with pytest.raises(ValueError):
            module.admit(tmp_path, submission_sha, protocol_sha, destination)
        assert not destination.exists()
    else:
        result = module.admit(tmp_path, submission_sha, protocol_sha, destination)
        assert result["status"] == "private_qfo_phylogeny_deployment_admitted_unscored"
        assert result["recovered_cpm_inference_authorized"] is True
        assert all(result[key] is False for key in ("accuracy_evaluated", "scoring_admitted", "controlled_timing", "publication_ready"))
        assert len(calls) == 2 and len(result["native_comparison"]) == 6
        assert result == json.loads(destination.read_bytes())
        with pytest.raises(FileExistsError):
            module.admit(tmp_path, submission_sha, protocol_sha, destination)


def test_output_symlink_refused(tmp_path):
    path = tmp_path / "admission.json"
    path.symlink_to(tmp_path / "missing")
    with pytest.raises(ValueError, match="direct absolute"):
        module.admit(tmp_path, "unused", "unused", path)


def test_cli_imports_no_scientific_package():
    code = ("import sys; import benchmark_tools.admit_qfo_private_phylogeny_control; "
            "assert not any(n == 'orthohmm' or n.startswith('orthohmm.') for n in sys.modules)")
    subprocess.run([sys.executable, "-S", "-B", "-c", code], check=True)


def test_script_pins_readonly_allocation():
    path = Path(module.__file__).parent / "results/qfo_private_phylogeny_admission_20261001.sh"
    subprocess.run(["bash", "-n", str(path)], check=True)
    text = path.read_text()
    for token in ("--cpus-per-task=2", "--mem=64G", "--time=04:00:00", "--nodelist=bizon", "--no-requeue",
                  "--submission-sha256", "--protocol-sha256"):
        assert token in text
    assert "--dependency" not in text and "--array" not in text and "dgx" not in text
