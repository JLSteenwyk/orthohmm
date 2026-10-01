import copy
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import admit_private_helper_cpm_phylogeny as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_admit_qfo_parameter_phylogeny import execution_fixture


def accounting(state="COMPLETED", exit_code="0:0", cpus="32", memory="192G", node="bizon", step="COMPLETED"):
    return ("JobID|JobIDRaw|State|ExitCode|Elapsed|AllocCPUS|ReqMem|NodeList\n"
        f"22390|22390|{state}|{exit_code}|01:00:00|{cpus}|{memory}|{node}\n"
        f"22390.batch|22390.batch|{step}|0:0|01:00:00|32||bizon\n")


@pytest.mark.parametrize("values", [{}, {"state": "RUNNING"}, {"state": "PENDING"}, {"state": "FAILED"},
    {"state": "TIMEOUT"}, {"state": "CANCELLED"}, {"exit_code": "1:0"}, {"cpus": "2"},
    {"memory": "64G"}, {"node": "dgx"}, {"step": "RUNNING"}, {"step": "FAILED"}])
def test_native_completion_gate(values):
    if values:
        with pytest.raises(ValueError):
            module.completed(accounting(**values))
    else:
        assert module.completed(accounting())["JobIDRaw"] == "22390"


@pytest.mark.parametrize("text", ["", accounting().split("22390.batch")[0],
    accounting() + accounting().splitlines()[1] + "\n", accounting() + accounting().splitlines()[2] + "\n"])
def test_missing_or_ambiguous_native_accounting(text):
    with pytest.raises(ValueError):
        module.completed(text)


def test_active_native_precedes_data_reads(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "accounting", lambda: accounting(state="RUNNING"))
    monkeypatch.setattr(module, "producer", lambda *args: pytest.fail("Read live native output"))
    with pytest.raises(ValueError, match="completed"):
        module.admit(tmp_path, "unused", "unused", tmp_path / "admission.json")
    assert not (tmp_path / "admission.json").exists()


@pytest.mark.parametrize("problem", [None, "post_extra", "post_missing", "publication", "numeric_flag", "command", "fresh", "status", "provenance"])
def test_exact_unscored_execution(problem):
    status, preflight, postflight, expected, manifest, cell = execution_fixture()
    command, fresh = ["historical-controller", "-B", "candidate-admitter.py"], {"path": "fresh.json", "sha256": "pin", "bytes": 1}
    postflight.update(publication_ready=False, admission_command=copy.deepcopy(command), fresh_candidate_admission=copy.deepcopy(fresh))
    if problem == "post_extra":
        postflight["error"] = "mixed failure/success"
    elif problem == "post_missing":
        del postflight["publication_ready"]
    elif problem == "publication":
        postflight["publication_ready"] = True
    elif problem == "numeric_flag":
        postflight["publication_ready"] = 0
    elif problem == "command":
        postflight["admission_command"][0] = "native-private-python"
    elif problem == "fresh":
        postflight["fresh_candidate_admission"]["sha256"] = "other"
    elif problem == "status":
        status["failed_methods"] = [cell["label"]]
    elif problem == "provenance":
        preflight["job_id"] = "other"
    if problem:
        with pytest.raises(ValueError):
            module.execution_identity(preflight, postflight, status, expected, manifest, cell, command, fresh)
    else:
        module.execution_identity(preflight, postflight, status, expected, manifest, cell, command, fresh)


def cache_fixture():
    return dict(candidate_families=346866, reconciled_families=100, bypassed_families=346766,
        checkpoint_hits=90, remapped_checkpoint_hits=5, species_tree_families=20, species_tree_checkpoint_hit=False)


@pytest.mark.parametrize("problem", [None, "candidate", "total", "bool", "negative", "hits", "remapped", "species_count", "species_hit"])
def test_actual_cache_accounting(problem):
    summary = cache_fixture()
    changes = {"candidate": ("candidate_families", 351739), "total": ("bypassed_families", 0),
        "bool": ("checkpoint_hits", True), "negative": ("checkpoint_hits", -1), "hits": ("checkpoint_hits", 101),
        "remapped": ("remapped_checkpoint_hits", 91), "species_count": ("species_tree_families", 346867),
        "species_hit": ("species_tree_checkpoint_hit", 0)}
    if problem:
        key, value = changes[problem]
        summary[key] = value
        with pytest.raises(ValueError):
            module.cache_accounting(summary)
    else:
        assert module.cache_accounting(summary) == summary


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return record(path)


@pytest.mark.parametrize("problem", [None, "returncode", "command", "request", "cwd", "executable", "bytecode", "prefix", "enforcement", "missing_pipeline", "file"])
def test_private_lookup_identity(tmp_path, problem):
    from benchmark_tools.inspect_native_python_lookup import PROBE, scientific_origins
    launcher, output = tmp_path / "launcher", tmp_path / "output"
    python = tmp_path / "private-python"
    python.write_text("private executable")
    for relative in ("orthohmm/phylogeny_pipeline.py", "benchmark_tools/replay_phylogeny.py"):
        path = launcher / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("frozen source")
    requested = ["orthohmm.phylogeny_pipeline", "benchmark_tools.replay_phylogeny"]
    report = dict(modules={name: str(launcher / relative) for name, relative in zip(requested,
        ("orthohmm/phylogeny_pipeline.py", "benchmark_tools/replay_phylogeny.py"))}, requested=requested,
        mapped_files=[], executable=str(python), cwd=str(launcher), dont_write_bytecode=True,
        pycache_prefix=str(output / "bytecode_cache"))
    paths = sorted({*report["modules"].values(), str(python)})
    lookup = dict(report=report, origin=scientific_origins(report, "orthohmm", launcher / "orthohmm"),
        continuous_enforcement=False, checked_records=[record(path) for path in paths])
    process = dict(returncode=0, command=[str(python), "-B", "-c", PROBE, json.dumps(requested)])
    changes = {"returncode": (process, "returncode", 1), "command": (process, "command", []),
        "request": (report, "requested", []), "cwd": (report, "cwd", "wrong"),
        "executable": (report, "executable", "wrong"), "bytecode": (report, "dont_write_bytecode", False),
        "prefix": (report, "pycache_prefix", "wrong"), "enforcement": (lookup, "continuous_enforcement", True)}
    if problem in changes:
        target, key, value = changes[problem]
        target[key] = value
    elif problem == "missing_pipeline":
        report["modules"].pop(requested[0])
    process["stdout"] = json.dumps(report)
    write_json(output / "lookup.json", lookup)
    write_json(output / "lookup_process.json", process)
    if problem == "file":
        python.write_text("changed executable")
    if problem:
        with pytest.raises(ValueError):
            module.lookup_identity(output, [str(python)], launcher)
    else:
        assert len(module.lookup_identity(output, [str(python)], launcher)) == 5


def admission_fixture(root, monkeypatch, problem):
    from benchmark_tools import run_simulation_methods as execution
    from benchmark_tools import validate_simulation_outputs as process_module
    from benchmark_tools import validate_factorial_native as native_module
    from benchmark_tools import admit_qfo_corrected_factorial_cell as ownership_module
    from benchmark_tools import admit_qfo_factorial_cell as pairs_module
    from benchmark_tools.verify_qfo_replay_launcher import LAUNCHER_COMMIT

    executor = root / module.EXECUTOR
    source = executor / "benchmark_tools/run_private_helper_cpm_phylogeny.py"
    source.parent.mkdir(parents=True)
    source.write_bytes((Path(module.__file__).parent / source.name).read_bytes())
    producer_protocol = root / module.PRODUCER_PROTOCOL
    producer_protocol.parent.mkdir(parents=True)
    producer_protocol.write_text("frozen producer protocol")
    executor_protocol = executor / module.PRODUCER_PROTOCOL
    executor_protocol.parent.mkdir(parents=True)
    executor_protocol.write_bytes(producer_protocol.read_bytes())
    monkeypatch.setattr(module, "PRODUCER_PROTOCOL_SHA", record(producer_protocol)["sha256"])
    protocol = root / module.PROTOCOL
    protocol.write_text("prospective admission protocol")
    launcher = root / "launcher"
    equivalents = []
    for name in ("replay_phylogeny.py", "orthobench_stage_diagnostics.py"):
        executed, prepared = launcher / "benchmark_tools" / name, root / "prepared/benchmark_tools" / name
        executed.parent.mkdir(parents=True, exist_ok=True)
        prepared.parent.mkdir(parents=True, exist_ok=True)
        executed.write_text("frozen launcher")
        prepared.write_bytes(executed.read_bytes())
        equivalents.append(dict(prepared=record(prepared), executed=record(executed)))
    extra = root / "evidence"
    extra.write_text("bound evidence")
    candidate = {}
    for key in ("seed_partition", "candidate_partition", "membership_constraints"):
        path = root / (key + ".txt")
        path.write_text(key)
        candidate[key] = record(path)
    arm = dict(seed_partition=candidate["seed_partition"], partition=candidate["candidate_partition"],
               constraints=candidate["membership_constraints"])
    candidate_report = dict(candidate_arm=candidate, preparation={"sha256": "preparation-pin"}, protocol={"sha256": "candidate-protocol-pin"})
    baseline = dict(launcher=str(launcher), manifest=dict(input_fastas=[], candidate_arms={"p1_c1": {"seed_partition": "old_seed"}}), environment={})
    verified = dict(baseline=baseline, candidates=dict(arm=arm, admission=candidate_report,
        admission_record=record(extra), admission_executor=str(root / "candidate-executor")),
        private_control=dict(admission=record(extra)), checked_records=[record(extra)])
    output_relative = "benchmarks/results/qfo_cpm_private_phylogeny_v1/cpm_high"
    output = root / output_relative
    argv = [str(root / "private-python"), str(launcher / "benchmark_tools/replay_phylogeny.py")]
    cell, planned = {"label": "candidate_cpm_high", "argv": argv}, {"label": "candidate_cpm_high", "argv": [sys.executable, "planned.py"]}
    calls = []
    def verify(*args):
        calls.append(args)
        return {} if problem == "changed" else verified
    fake = SimpleNamespace(OUTPUT=output_relative, PYTHON=sys.executable,
        native_command=lambda *args: (cell, planned, argv, equivalents), verify_sources=verify)
    monkeypatch.setattr(module, "producer", lambda *args: (fake, verified))
    monkeypatch.setattr(module, "PRODUCER_SHA", "wrong" if problem == "source" else record(source)["sha256"])
    monkeypatch.setattr(module, "accounting", lambda: accounting())
    monkeypatch.setattr(module.subprocess, "check_output", lambda *args, **kwargs: "wrong" if problem == "revision" else module.COMMIT)
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: None)
    monkeypatch.setattr(execution, "execution_environment", lambda *args: ({}, {}))
    expected = dict(source=record(source), protocol=record(producer_protocol), helpers=[record(source)], verified=verified,
        cell=cell, planned_cell=planned, executed_argv=argv, launcher_source_equivalence=equivalents,
        resolved_tools={}, cwd=str(launcher), job_id=module.JOB,
        seed_handoff="explicit_helper_runtime_seed_amendment", native_handoff="explicit_admitted_private_phylogeny_deployment",
        scope="Unscored fixed recovered high-CPM arm; inferred phylogeny with validated raw-tree checkpoint reuse; incremental shared-host execution")
    preflight = copy.deepcopy(expected)
    if problem == "helpers":
        preflight["helpers"] = []
    elif problem == "protocol_path":
        preflight["protocol"] = record(executor_protocol)
    fresh = write_json(output / "fresh_candidate_admission.json", {} if problem == "fresh" else candidate_report)
    command = [sys.executable, "-B", str(root / "candidate-executor/benchmark_tools/admit_helper_cpm_candidates.py"),
        "--root", str(root), "--preparation-sha256", "preparation-pin", "--protocol-sha256", "candidate-protocol-pin", "--output", fresh["path"]]
    postflight = dict(status="complete_pending_native_validation", cell=cell, accuracy_evaluated=False,
        native_outputs_validated=False, publication_ready=False, admission_command=command, fresh_candidate_admission=fresh)
    status = dict(provenance=preflight, verified_inputs=dict(status="ready", inputs=[]), dataset=cell["label"],
        methods={cell["label"]: {}}, status="finished_pending_native_validation", failed_methods=[],
        accuracy_evaluated=False, native_outputs_validated=False)
    if problem == "execution":
        status["status"] = "running"
    for relative, value in (("preflight.json", preflight), ("postflight.json", postflight), ("execution/status.json", status)):
        write_json(output / relative, value)
    (output / "admission.log").write_text("")
    native_dir = output / "output/orthohmm_phylogeny"
    write_json(native_dir / "provenance_manifest.json", dict(results={**cache_fixture(), "ortholog_pairs": 1}))
    write_json(output / "metrics.json", {})
    native_dir.joinpath("species_tree.rooted.nwk").write_text("(a,b);\n")
    native_dir.joinpath("orthohmm_pairwise_orthologs.tsv").write_text("a\tb\n")
    monkeypatch.setattr(module, "lookup_identity", lambda *args: [])
    events = []
    monkeypatch.setattr(process_module, "verify_process", lambda *args: events.append("process") or set())
    def native(adapted, *args, **kwargs):
        events.append("native")
        assert adapted["candidate_arms"]["p1_c1"] == candidate
        assert kwargs["expected_revision"] == LAUNCHER_COMMIT
        return dict(native_manifest=record(native_dir / "provenance_manifest.json"), native_metrics=record(output / "metrics.json"), species_tree=record(native_dir / "species_tree.rooted.nwk"))
    monkeypatch.setattr(native_module, "validate_native_cell", native)
    monkeypatch.setattr(ownership_module, "gene_ownership", lambda *args: ({}, {}))
    monkeypatch.setattr(pairs_module, "check_pairs", lambda *args: 0 if problem == "pairs" else 1)
    submission = dict(status="private_recovered_qfo_high_cpm_phylogeny_submitted", job_id=module.JOB,
        executor=str(executor), executor_commit=module.COMMIT, executor_clean=True, source_records=[record(source)],
        source_preflight=record(extra), private_admission=record(extra), private_readback=record(extra), candidate_admission=record(extra),
        native_completion_observed=False, native_outputs_validated=False, accuracy_evaluated=False, controlled_timing=False, publication_ready=False)
    submission_pin = write_json(root / module.SUBMISSION, submission)
    return submission_pin["sha256"], record(protocol)["sha256"], calls, events


@pytest.mark.parametrize("problem", [None, "revision", "source", "helpers", "protocol_path", "fresh", "execution", "pairs", "changed"])
def test_complete_admission_contract(tmp_path, monkeypatch, problem):
    submission_sha, protocol_sha, calls, events = admission_fixture(tmp_path, monkeypatch, problem)
    destination = tmp_path / "admission.json"
    if problem:
        with pytest.raises(ValueError):
            module.admit(tmp_path, submission_sha, protocol_sha, destination)
        assert not destination.exists()
    else:
        result = module.admit(tmp_path, submission_sha, protocol_sha, destination)
        assert result["status"] == "private_recovered_cpm_native_pairs_verified_unscored"
        assert result["arm"] == "cpm_high" and result["index"] == 1 and result["native_pair_count"] == 1
        assert all(result[k] is False for k in ("accuracy_evaluated", "scoring_admitted", "controlled_timing", "publication_ready"))
        assert len(calls) == 1 and events == ["process", "native", "process"]
        assert result == json.loads(destination.read_bytes())
        with pytest.raises(FileExistsError):
            module.admit(tmp_path, submission_sha, protocol_sha, destination)


def test_symlink_report_refused(tmp_path):
    destination = tmp_path / "admission.json"
    destination.symlink_to(tmp_path / "missing")
    with pytest.raises(ValueError, match="direct absolute"):
        module.admit(tmp_path, "unused", "unused", destination)


def test_cli_does_not_import_scientific_code():
    code = ("import sys; import benchmark_tools.admit_private_helper_cpm_phylogeny; "
            "assert not any(n == 'orthohmm' or n.startswith('orthohmm.') for n in sys.modules)")
    subprocess.run([sys.executable, "-S", "-B", "-c", code], check=True)


def test_readonly_script_contract():
    path = Path(module.__file__).parent / "results/qfo_private_cpm_native_admission_20261001.sh"
    subprocess.run(["bash", "-n", str(path)], check=True)
    text = path.read_text()
    for token in ("--cpus-per-task=2", "--mem=64G", "--time=04:00:00", "--no-requeue", "--nodelist=bizon",
                  "--submission-sha256", "--protocol-sha256", "admit_private_helper_cpm_phylogeny.py"):
        assert token in text
    assert "--dependency" not in text and "--array" not in text and "dgx" not in text
