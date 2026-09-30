import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import helper_recovered_cpm_seed_evidence as module
from benchmark_tools import prepare_recovered_cpm_candidates as candidates
from benchmark_tools.prepare_ob_candidate_neighborhood import record


ACCOUNTING = """JobID|JobIDRaw|State|ExitCode|Elapsed|AllocCPUS|ReqMem|NodeList
22154|22154|COMPLETED|0:0|00:11:38|1|64G|bizon
22081_1|22081|FAILED|1:0|00:48:46|32|192G|bizon
22155|22155|FAILED|1:0|00:01:47|2|64G|bizon
22156|22156|CANCELLED by 1000|0:0|00:00:00|0|64G|None assigned
"""


@pytest.mark.parametrize("problem", [None, "readback_sha", "selected_sha", "protocol_sha", "report_sha",
    "admission_record", "source", "revision", "status", "readback_status", "seed_flag", "accuracy",
    "downstream", "publication", "native", "scope", "optimizer_runtime", "downstream_runtime", "stages",
    "seed", "coverage", "graph", "constructor", "record_count", "missing_provenance", "changed_seed",
    "changed_transitive", "missing_binding", "binding_revision", "binding_blob", "native_pending",
    "failure_relabelled", "old_candidate_pending", "scheduler_evidence"])
def test_explicit_seed_handoff(tmp_path, monkeypatch, problem):
    def file(relative, content="fixture"):
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return record(path)
    source = file("benchmark_tools/admit_helper_cpm_recovery.py")
    monkeypatch.setattr(module, "SOURCE_SHA", "changed" if problem == "source" else source["sha256"])
    protocol = file(module.PROTOCOL)
    relative = ("benchmark_tools/admit_helper_cpm_recovery.py", "tests/unit/test_admit_helper_cpm_recovery.py",
        "benchmark_tools/results/QFO_CPM_HELPER_RECOVERY_ADMISSION_PROTOCOL_20260930.md",
        "benchmark_tools/probe_cpm_partition_parser.py", "benchmark_tools/probe_cpm_private_runtime.py")
    bindings = [dict(file(path), git_revision=module.REVISION) for path in relative]
    parent = file("benchmarks/results/qfo_cpm_checkpoint_recovery_v1/status.json")
    original = "benchmarks/results/qfo_parameter_cpm_replay_v3/cpm_high/replay/"
    recovered = "benchmarks/results/qfo_cpm_checkpoint_recovery_v1/"
    outputs = [file(path, "a b\n") for path in (original + "orthogroups_multipass.txt",
        original + "orthogroups_multipass_refined.txt", recovered + "orthogroups_profiles.txt",
        recovered + "orthogroups_profiles_refined.txt")]
    stages = [dict(label=label, origin=origin, output=output) for label, origin, output in zip(
        ("multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"),
        ("reused", "reused", "recovered", "recovered"), outputs)]
    coverage = [dict(label=stage["label"], genes=984137, groups=groups) for stage, groups in
                zip(stages, (316603, 393142, 314274, 390845))]
    scheduler = module.scheduler_rows(ACCOUNTING)
    transitive = file("transitive.txt")
    records = [source, parent, transitive, *outputs]
    report = dict(status="cpm_helper_runtime_recovered_seed_admitted_unscored", source=source,
        seed_admitted=True, native_attempts=0, accuracy_evaluated=False, downstream_admitted=False,
        publication_ready=False, stages=stages, seed_partition=outputs[-1], source_report=parent,
        coverage=[dict(row, output=stage["output"]) for row, stage in zip(coverage, stages)],
        saved_graph={"fixture": "graph"}, constructor_bytes_sha256="constructor",
        checked_records=records, scheduler=scheduler["22154"], original_failure=scheduler["22081_1"],
        historical_failed_admission=scheduler["22155"], refinement_runtime_amendment=dict(
            scope="independent refinement only", original_optimizer_runtime_unchanged=True,
            downstream_runtime_authorized=False))
    readback = dict(status="helper_runtime_recovered_seed_admission_independently_read_back_unscored",
        source_revision=module.REVISION, source_git_bindings=bindings, seed_admitted=True, native_attempts=0,
        accuracy_evaluated=False, downstream_admitted=False, publication_ready=False, coverage=coverage,
        seed_partition=outputs[-1], saved_graph=report["saved_graph"], constructor_bytes_sha256="constructor",
        checked_records_reverified=len(records), scheduler=report["scheduler"],
        original_failure=report["original_failure"], historical_failed_admission=report["historical_failed_admission"])
    if problem in ("status", "seed_flag", "accuracy", "downstream", "publication", "native", "stages", "seed", "missing_provenance"):
        key, value = {"status": ("status", "running"), "seed_flag": ("seed_admitted", False),
            "accuracy": ("accuracy_evaluated", True), "downstream": ("downstream_admitted", True),
            "publication": ("publication_ready", True), "native": ("native_attempts", 1),
            "stages": ("stages", list(reversed(stages))), "seed": ("seed_partition", outputs[2]),
            "missing_provenance": ("checked_records", [])}[problem]
        report[key] = value
    elif problem in ("scope", "optimizer_runtime", "downstream_runtime"):
        key, value = {"scope": ("scope", "all science"), "optimizer_runtime": ("original_optimizer_runtime_unchanged", False),
                      "downstream_runtime": ("downstream_runtime_authorized", True)}[problem]
        report["refinement_runtime_amendment"][key] = value
    elif problem == "scheduler_evidence":
        report["scheduler"] = dict(report["scheduler"], Elapsed="different")
    elif problem == "revision":
        readback["source_revision"] = "different"
    elif problem == "readback_status":
        readback["status"] = "running"
    elif problem == "coverage":
        readback["coverage"] = []
    elif problem == "graph":
        readback["saved_graph"] = {}
    elif problem == "constructor":
        readback["constructor_bytes_sha256"] = "different"
    elif problem == "record_count":
        readback["checked_records_reverified"] += 1
    elif problem == "missing_binding":
        readback["source_git_bindings"] = bindings[:-1]
    elif problem == "binding_revision":
        bindings[0]["git_revision"] = "different"
    report_record = file(module.ADMISSION, json.dumps(report))
    monkeypatch.setattr(module, "ADMISSION_SHA", "changed" if problem == "report_sha" else report_record["sha256"])
    readback["admission"] = dict(report_record, sha256="changed") if problem == "admission_record" else report_record
    readback_record = file(module.READBACK, json.dumps(readback))
    monkeypatch.setattr(module, "READBACK_SHA", "changed" if problem == "readback_sha" else readback_record["sha256"])
    if problem == "changed_seed":
        Path(outputs[-1]["path"]).write_text("changed")
    elif problem == "changed_transitive":
        Path(transitive["path"]).write_text("changed")
    accounting = ACCOUNTING
    if problem == "native_pending":
        accounting = accounting.replace("22154|22154|COMPLETED", "22154|22154|PENDING")
    elif problem == "failure_relabelled":
        accounting = accounting.replace("22155|22155|FAILED", "22155|22155|COMPLETED")
    elif problem == "old_candidate_pending":
        accounting = accounting.replace("CANCELLED by 1000", "PENDING")
    def check_output(command, **kwargs):
        if command[0] == "sacct":
            return accounting
        assert command[:4] == ["git", "-C", str(tmp_path), "show"]
        return b"changed" if problem == "binding_blob" else b"fixture"
    monkeypatch.setattr(module.subprocess, "check_output", check_output)
    selection = (tmp_path, Path(readback_record["path"]),
        "different" if problem == "selected_sha" else readback_record["sha256"],
        "different" if problem == "protocol_sha" else protocol["sha256"])
    if problem:
        with pytest.raises(ValueError):
            module.evidence(*selection)
    else:
        result = module.evidence(*selection)
        assert result["seed_partition"] == outputs[-1]
        assert result["handoff"] == "explicit_helper_runtime_seed_amendment"
        assert result["protocol"] == protocol
        assert result["accuracy_evaluated"] is False
        assert result["scheduler"]["22155"]["State"] == "FAILED"
        assert result["scheduler"]["22156"]["State"] == "CANCELLED by 1000"


@pytest.mark.parametrize("change", ["duplicate", "missing", "memory", "cpu", "exit", "host", "candidate_executed"])
def test_scheduler_contract_rejects(change):
    accounting = ACCOUNTING
    if change == "duplicate":
        accounting += ACCOUNTING.splitlines()[1] + "\n"
    elif change == "missing":
        accounting = "\n".join(row for row in accounting.splitlines() if not row.startswith("22154|"))
    elif change == "memory":
        accounting = accounting.replace("|1|64G|bizon", "|1|32G|bizon")
    elif change == "cpu":
        accounting = accounting.replace("|1|64G|bizon", "|2|64G|bizon")
    elif change == "exit":
        accounting = accounting.replace("22154|22154|COMPLETED|0:0", "22154|22154|COMPLETED|1:0")
    elif change == "host":
        accounting = accounting.replace("bizon", "elsewhere")
    elif change == "candidate_executed":
        accounting = accounting.replace("00:00:00|0|64G", "00:01:00|2|64G")
    with pytest.raises(ValueError):
        module.scheduler_rows(accounting)


def test_partial_cli_selection_is_argparse_error(tmp_path):
    done = subprocess.run([sys.executable, "-I", str(Path(candidates.__file__).resolve()),
        "--root", str(tmp_path), "--helper-recovery-readback", "fixture"], capture_output=True, text=True, cwd=tmp_path)
    assert done.returncode == 2
    assert "required together" in done.stderr
    assert not (tmp_path / "benchmarks").exists()


def test_historical_gate_remains_unchanged():
    source = Path(module.__file__).with_name("recovered_cpm_seed_evidence.py")
    assert record(source)["sha256"] == "7e47f010ffd6340f17074e66e0f29836d48f16db543c94d8f20b88279b6b1046"


def test_launch_script_has_explicit_amendment_and_original_runtime():
    root = Path(module.__file__).resolve().parent.parent
    script = root / "benchmark_tools/results/qfo_cpm_helper_candidates_20260930.sh"
    done = subprocess.run(["bash", "-n", str(script)], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    content = script.read_text()
    assert "--helper-recovery-readback-sha256 " + module.READBACK_SHA in content
    assert "--helper-candidate-protocol-sha256 " + record(root / module.PROTOCOL)["sha256"] in content
    assert "/home/bizon/anaconda3/bin/python -B" in content
    for flag in ("--nodelist=bizon", "--ntasks=1", "--cpus-per-task=2", "--mem=64G", "--time=04:00:00", "--no-requeue"):
        assert "#SBATCH " + flag in content
    assert "--dependency" not in content
    assert "dgx" not in content.lower()
