from copy import deepcopy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from benchmark_tools import assess_three_kingdoms_matched_sonic as module


def synthetic_reference_groups():
    # Match only the admission dimensions, never BUSCO identities or families.
    sizes = [29] + [7] * 10 + [6] * 4 + [4] * 2 + [8] * 238
    return [[f"synthetic_family_{family}_gene_{gene}" for gene in range(size)]
            for family, size in enumerate(sizes)]


def test_synthetic_reference_dimensions_and_distinct_identities():
    groups = synthetic_reference_groups()
    genes = [gene for group in groups for gene in group]
    assert len(groups) == 255
    assert len(genes) == len(set(genes)) == 2035
    assert sum(len(group) * (len(group) - 1) // 2 for group in groups) == 7352
    assert all(gene.startswith("synthetic_family_") for gene in genes)


@pytest.fixture
def case():
    plan = {"output_root": "/run", "native_argv": ["sonic", "-t", "32"],
            "inputs": [{"path": "/input/a.fasta", "bytes": 10, "sha256": "abc"}]}
    execution = {"status": "process_succeeded_pending_native_admission", "exit_code": 0,
                 "accuracy_admitted": False, "job_id": "21795", "node": "bizon",
                 "plan": {"sha256": "plan"}, "source": {"sha256": "runner"},
                 "native_argv": plan["native_argv"], "started_epoch": 1, "finished_epoch": 2,
                 "runtime_before": {"tree": "ok"}, "runtime_after": {"tree": "ok"},
                 "copied_inputs": [{"path": "/run/input/a.fasta", "bytes": 10, "sha256": "abc"}],
                 "native_log": {"path": "/run/native.log", "bytes": 1},
                 "timing": {"path": "/run/time.txt", "bytes": 1}}
    scheduler = {"JobIDRaw": "21795", "State": "COMPLETED", "ExitCode": "0:0",
                 "NodeList": "bizon", "AllocCPUS": "32", "ReqMem": "192G"}
    return plan, execution, scheduler, deepcopy(execution["plan"]), deepcopy(execution["source"])


def test_valid_execution(case):
    module.validate_execution(*case)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
    ("JobIDRaw", "21796"), ("AllocCPUS", "16"), ("ReqMem", "64G"), ("NodeList", "other")])
def test_wrong_scheduler_rejected(case, key, value):
    case[2][key] = value
    with pytest.raises(ValueError):
        module.validate_execution(*case)


@pytest.mark.parametrize("key,value", [("status", "running"), ("exit_code", 1),
    ("accuracy_admitted", True), ("job_id", "other"), ("node", "other"),
    ("plan", {}), ("source", {}), ("native_argv", []), ("started_epoch", 3),
    ("runtime_after", {}), ("copied_inputs", []), ("timing", {"path": "/bad", "bytes": 1})])
def test_wrong_execution_rejected(case, key, value):
    case[1][key] = value
    with pytest.raises(ValueError):
        module.validate_execution(*case)


def test_pending_job_cannot_create_assessment(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
                        "JobIDRaw|State|ExitCode\n21795|PENDING|0:0\n")
    destination = tmp_path / "new"
    with pytest.raises(ValueError, match="COMPLETED"):
        module.assess(tmp_path, destination)
    assert not destination.exists()


def test_existing_destination_preserved(tmp_path):
    with pytest.raises(FileExistsError):
        module.assess(tmp_path, tmp_path)


@pytest.fixture(params=["synthetic_dimensions", "retained_busco_reference"])
def integrated(tmp_path, monkeypatch, request):
    """Real scoring on synthetic and, when available, pinned retained references."""
    repo = tmp_path / "repo"
    source = Path(__file__).resolve().parents[2]
    for name in ("benchmark_tools/normalize_three_kingdoms_orthogroups.py",
                 "three_kingdoms/score_against_busco.py"):
        target = repo / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source / name, target)
    reference = repo / "three_kingdoms/busco/reference_orthogroups.txt"
    reference.parent.mkdir(parents=True)
    if request.param == "retained_busco_reference":
        retained = source / "three_kingdoms/busco/reference_orthogroups.txt"
        if not retained.is_file():
            pytest.skip("Retained BUSCO data unavailable; synthetic integration is separate")
        historical = json.loads((source / "benchmark_tools/results/"
            "three_kingdoms_sonic_matched_commands_20260918.json").read_text())
        pins = [item for item in historical["checked_records"]
                if item["path"].endswith("/three_kingdoms/busco/reference_orthogroups.txt")]
        assert len(pins) == 1
        module.check({**pins[0], "path": str(retained)})
        shutil.copyfile(retained, reference)
    else:
        reference.write_text("".join(" ".join(group) + "\n"
                                    for group in synthetic_reference_groups()))
    groups = [line.split() for line in reference.read_text().splitlines() if line.strip()]
    members = [[], []]
    rows = ["group_id\tgroup_size\tsp_in_grp\tseed_ortholog_cnt\ta.fasta\tb.fasta"]
    for index, genes in enumerate(groups, 1):
        cells = [genes[::2], genes[1::2]]
        for bucket, cell in zip(members, cells):
            bucket.extend(cell)
        rows.append("\t".join([str(index), str(len(genes)), str(sum(bool(c) for c in cells)),
                               str(len(genes)), *[",".join(c) or "*" for c in cells]]))
    root = repo / "native"
    inputs = root / "input"
    outputs = root / "output"
    inputs.mkdir(parents=True)
    outputs.mkdir()
    copies = []
    snapshot = []
    for index, (name, genes) in enumerate(zip(("a", "b"), members), 1):
        path = inputs / (name + ".fasta")
        path.write_text("".join(f">{gene}\nACDE\n" for gene in genes))
        item = module.record(path)
        copies.append(item)
        snapshot.append(f"{index}\t{path.name}\t{item['sha256']}\t{len(genes)}\t{4*len(genes)}")
    (outputs / "snapshot.tsv").write_text("\n".join(snapshot) + "\n")
    (outputs / "ortholog_groups.tsv").write_text("\n".join(rows) + "\n")
    (root / "native.log").write_text("synthetic native log\n")
    (root / "time.txt").write_text("synthetic timing, not a measurement\n")
    executor = repo / "benchmarks/work/publication_three_kingdoms_sonic_matched_v1"
    runner = executor / "benchmark_tools/run_three_kingdoms_matched_sonic.py"
    runner.parent.mkdir(parents=True)
    runner.write_text("# Synthetic executor identity fixture.\n")
    runtime_path = repo / "runtime.json"
    runtime_path.write_text(json.dumps({"python": {"path": sys.executable}}))
    plan = {"output_root": str(root), "native_argv": ["synthetic-inference-never-executed"],
            "inputs": copies, "runtime": module.record(runtime_path),
            "checked_records": [module.record(repo / name) for name in
                ("benchmark_tools/normalize_three_kingdoms_orthogroups.py",
                 "three_kingdoms/score_against_busco.py",
                 "three_kingdoms/busco/reference_orthogroups.txt")]}
    plan_path = repo / "benchmark_tools/results/three_kingdoms_sonic_matched_commands_20260918.json"
    plan_path.parent.mkdir()
    plan_path.write_text(json.dumps(plan))
    env = {"PATH": os.defpath, "PYTHONDONTWRITEBYTECODE": "1", "PYTHONNOUSERSITE": "1"}
    runtime = {"synthetic": "host inventory is mocked"}
    execution = {"status": "process_succeeded_pending_native_admission", "exit_code": 0,
                 "accuracy_admitted": False, "job_id": str(module.JOB), "node": "bizon",
                 "plan": module.record(plan_path), "source": module.record(runner),
                 "native_argv": plan["native_argv"], "started_epoch": 1, "finished_epoch": 2,
                 "runtime_before": runtime, "runtime_after": runtime, "environment": env,
                 "copied_inputs": copies, "native_log": module.record(root / "native.log"),
                 "timing": module.record(root / "time.txt"),
                 "outputs": [module.record(p) for p in sorted(outputs.iterdir())]}
    (root / "execution.json").write_text(json.dumps(execution))
    monkeypatch.setattr(module, "PLAN_SHA", module.record(plan_path)["sha256"])
    monkeypatch.setattr(module, "verify", lambda *args: (plan, root, env, runtime))
    original_run = subprocess.run

    def output(args, **kwargs):
        if args[0] == "sacct":
            return ("JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n"
                    f"{module.JOB}|COMPLETED|0:0|00:00:01|bizon|32|192G\n")
        assert args[0] == "git"
        return module.EXECUTOR + "\n"

    def run(args, **kwargs):
        if args[0] == "git":
            return subprocess.CompletedProcess(args, 0)
        return original_run(args, **kwargs)

    monkeypatch.setattr(module.subprocess, "check_output", output)
    monkeypatch.setattr(module.subprocess, "run", run)
    return repo, tmp_path / "assessment", run


def test_real_conversion_scoring_and_independent_audit(integrated):
    repo, destination, _ = integrated
    result = module.assess(repo, destination)
    assert result["status"] == "matched_three_kingdoms_sonic_score_verified"
    assert result["accuracy_admitted"] is True
    assert result["counts"]["true_positive_gene_pairs"] == 7352
    assert result["counts"]["false_positive_gene_pairs"] == 0
    assert result["counts"]["false_negative_gene_pairs"] == 0
    assert result["counts"]["f_score"] == 1.0
    assert result["native_validation"]["grouped_genes"] == 2035
    assert json.loads((destination / "report.json").read_text()) == result


def test_valid_but_imperfect_predictions_are_admitted(integrated):
    repo, destination, _ = integrated
    table = repo / "native/output/ortholog_groups.tsv"
    rows = table.read_text().splitlines()
    fields = rows[1].split("\t")
    first_species = fields[4].split(",")
    separated = first_species.pop()
    fields[1] = fields[3] = str(int(fields[1]) - 1)
    fields[4] = ",".join(first_species)
    rows[1] = "\t".join(fields)
    rows.append(f"extra\t1\t1\t1\t{separated}\t*")
    table.write_text("\n".join(rows) + "\n")
    execution_path = repo / "native/execution.json"
    execution = json.loads(execution_path.read_text())
    execution["outputs"] = [module.record(p) for p in sorted(table.parent.iterdir())]
    execution_path.write_text(json.dumps(execution))
    result = module.assess(repo, destination)
    assert result["accuracy_admitted"] is True
    assert 0 < result["counts"]["f_score"] < 1
    lost_pairs = int(fields[1])
    assert result["counts"]["false_negative_gene_pairs"] == lost_pairs
    assert result["counts"]["true_positive_gene_pairs"] == 7352 - lost_pairs
    assert result["counts"]["false_positive_gene_pairs"] == 0
    assert result["counts"]["f_score"] == pytest.approx(2 * (7352 - lost_pairs) / (14704 - lost_pairs))


def test_changed_reference_rejected_before_assessment(integrated):
    repo, destination, _ = integrated
    reference = repo / "three_kingdoms/busco/reference_orthogroups.txt"
    reference.write_text(reference.read_text() + "\n")
    with pytest.raises(ValueError, match="Frozen input/source identity changed"):
        module.assess(repo, destination)
    assert not destination.exists()


def test_repinned_reference_universe_change_still_rejected(integrated, monkeypatch):
    repo, destination, _ = integrated
    reference = repo / "three_kingdoms/busco/reference_orthogroups.txt"
    reference.write_text("\n".join(reference.read_text().splitlines()[1:]) + "\n")
    path = repo / "benchmark_tools/results/three_kingdoms_sonic_matched_commands_20260918.json"
    plan = json.loads(path.read_text())
    plan["checked_records"] = [module.record(reference) if item["path"] == str(reference)
                               else item for item in plan["checked_records"]]
    path.write_text(json.dumps(plan))
    monkeypatch.setattr(module, "PLAN_SHA", module.record(path)["sha256"])
    execution_path = repo / "native/execution.json"
    execution = json.loads(execution_path.read_text())
    execution["plan"] = module.record(path)
    execution_path.write_text(json.dumps(execution))
    with pytest.raises(ValueError, match="Reference universe changed"):
        module.assess(repo, destination)
    report = json.loads((destination / "report.json").read_text())
    assert report["status"] == "assessment_failed" and report["accuracy_admitted"] is False
    assert (destination / "score.txt").is_file()


@pytest.mark.parametrize("stage", ["conversion", "scoring"])
def test_corrupt_pipeline_output_retains_failure(integrated, monkeypatch, stage):
    repo, destination, run = integrated

    def corrupt(args, **kwargs):
        completed = run(args, **kwargs)
        if stage == "conversion" and "sonicparanoid" in args:
            path = destination / "orthogroups.txt"
            path.write_text("\n".join(path.read_text().split()) + "\n")
        elif stage == "scoring" and "--predictions" in args:
            kwargs["stdout"].flush()
            path = destination / "score.txt"
            path.write_text(path.read_text().replace("7,352", "7,351"))
        return completed

    monkeypatch.setattr(module.subprocess, "run", corrupt)
    with pytest.raises(ValueError):
        module.assess(repo, destination)
    report = json.loads((destination / "report.json").read_text())
    assert report["status"] == "assessment_failed"
    assert report["accuracy_admitted"] is False
    assert "error" in report
