import json
from types import SimpleNamespace

import pytest

from benchmark_tools import run_qfo_cpm_assessment as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture(index=0):
    recheck = {"path": "/recheck", "bytes": 1, "sha256": "check"}
    stage = {"status": "cpm_native_pairs_prepared_unscored", "arm": module.ARMS[index],
        "index": index, "participant": "ohmm_qfo_parameter_" + module.ARMS[index],
        "semantics": "native phylogenetically inferred pairs", "job_id": "23000", "array_task_id": str(index),
        "accuracy_evaluated": False, "publication_ready": False, "written_pairs": 2,
        "total_pairs": 2, "retained_pairs": 2, "removed_mapping_pairs": 0,
        "pairs": {"bytes": 8, "sha256": "pairs"}, "filtered_pairs": {"bytes": 8, "sha256": "pairs"},
        "native_admission_recheck": recheck, "checked_records": [recheck]}
    return stage, {"JobIDRaw": "23000"}


@pytest.mark.parametrize("index", range(2))
def test_both_prespecified_arms(index):
    stage, scheduler = fixture(index)
    module.validate_stage(stage, index, scheduler)


@pytest.mark.parametrize("key,value", [("status", "other"), ("arm", "control"), ("index", 1),
    ("index", False), ("participant", "other"), ("semantics", "group cliques"), ("job_id", "999"),
    ("array_task_id", "1"), ("accuracy_evaluated", True), ("publication_ready", True),
    ("written_pairs", 3), ("total_pairs", 3), ("retained_pairs", True), ("removed_mapping_pairs", 1),
    ("filtered_pairs", {"bytes": 9, "sha256": "pairs"}), ("checked_records", [])])
def test_reject_changed_conversion(key, value):
    stage, scheduler = fixture()
    stage[key] = value
    with pytest.raises(ValueError):
        module.validate_stage(stage, 0, scheduler)


@pytest.mark.parametrize("index", [-1, 2, True, "0", 0.0])
def test_invalid_index(index):
    with pytest.raises(ValueError):
        module.validate_stage({}, index, {})
    with pytest.raises(ValueError):
        module.completed_conversion("", index)
    with pytest.raises(ValueError):
        module.prepare(None, index)


@pytest.mark.parametrize("state", ["RUNNING", "FAILED", "PENDING"])
def test_terminal_gate_precedes_artifact_read(tmp_path, monkeypatch, state):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        f"JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS\n22072_0|23000|{state}|0:0|bizon|2\n")
    monkeypatch.setattr(module, "record", lambda *a: pytest.fail("Read unfinished conversion"))
    with pytest.raises(ValueError, match="terminal"):
        module.prepare(tmp_path, 0)


@pytest.mark.parametrize("row", ["22072_0|23000|COMPLETED|1:0|bizon|2",
    "22072_0|23000|COMPLETED|0:0|bizon|8", "22072_0|23000|COMPLETED|0:0|other|2",
    "22072_0.batch|23000|COMPLETED|0:0|bizon|2", "22072_1|23000|COMPLETED|0:0|bizon|2"])
def test_wrong_scheduler_identity(row):
    with pytest.raises(ValueError, match="terminal"):
        module.completed_conversion("JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS\n" + row + "\n", 0)


@pytest.mark.parametrize("index", range(2))
def test_successful_raw_job_binding(index):
    row = module.completed_conversion(
        f"JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS\n22072_{index}|23000|COMPLETED|0:0|bizon|2\n", index)
    assert row["JobIDRaw"] == "23000"


@pytest.mark.parametrize("outcome", ["success", "exit_failure", "exception", "source_changed"])
def test_scoring_retains_evidence_without_admitting_accuracy(tmp_path, monkeypatch, outcome):
    source = tmp_path / "source.py"
    source.write_text("# frozen\n")
    results = tmp_path / "results"
    results.mkdir()
    (results / "score.json").write_text("{}")
    output = tmp_path / "execution"
    prepared = {"status": "prepared_unrun", "cwd": str(output), "results": str(results),
        "source": record(source), "verified_records": [], "environment_overrides": {"OMP_NUM_THREADS": "1"},
        "command": ["frozen-nextflow", "six-endpoints"], "accuracy_admitted": False}
    monkeypatch.setattr(module, "prepare", lambda *a: prepared.copy())
    for key, value in {"SLURM_JOB_ID": "23100", "SLURM_CPUS_PER_TASK": "8",
                       "SLURM_JOB_NODELIST": "bizon", "SLURM_ARRAY_TASK_ID": "0"}.items():
        monkeypatch.setenv(key, value)

    def execute(argv, **kwargs):
        assert argv == prepared["command"]
        assert kwargs["cwd"] == output
        assert kwargs["env"]["OMP_NUM_THREADS"] == "1"
        if outcome == "exception":
            raise RuntimeError("Interrupted")
        if outcome == "source_changed":
            source.write_text("# changed\n")
        return SimpleNamespace(returncode=1 if outcome == "exit_failure" else 0)

    monkeypatch.setattr(module.subprocess, "run", execute)
    if outcome == "success":
        module.run(tmp_path, 0)
    else:
        with pytest.raises((RuntimeError, ValueError)):
            module.run(tmp_path, 0)
    report = json.loads((output / "results.json").read_text())
    assert report["status"] == ("process_succeeded_pending_independent_admission" if outcome == "success" else "failed")
    assert report["accuracy_admitted"] is False
    assert report["array_task_id"] == "0"
    assert (output / "preflight.json").exists()
    with pytest.raises(FileExistsError):
        module.run(tmp_path, 0)


def test_allocation_gate_precedes_preparation(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    monkeypatch.setattr(module, "prepare", lambda *a: pytest.fail("Prepared without allocation"))
    with pytest.raises(ValueError, match="scheduled"):
        module.run(tmp_path, 0)


@pytest.mark.parametrize("change", [None, "context", "fastas", "source", "mapping", "counts",
                                   "output", "work", "results", "symlink", "payload"])
def test_prepare_binds_real_records_and_fresh_namespaces(tmp_path, monkeypatch, change):
    def write(relative, content):
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return record(path)

    source = write("benchmarks/work/publication_qfo_cpm_pairs_v2/benchmark_tools/prepare_qfo_cpm_pairs.py", "# converter\n")
    monkeypatch.setattr(module, "CONVERTER_SHA", source["sha256"])
    mapping = write("reference/mapping.json.gz", "mapping")
    pairs = write("pairs.tsv", "a\tb\n")
    recheck = write("native_admission_recheck.json", "{}")
    env = {"reference_files": [mapping], "environment_overrides": {"OMP_NUM_THREADS": "1"}}
    env_record = write("benchmark_tools/results/qfo_assessment_environment_20260917.json", json.dumps(env))
    monkeypatch.setattr(module, "ENV_SHA", env_record["sha256"])
    stage, _ = fixture()
    stage.update(source=source, context={"arm": "cpm_low", "resolution": .08}, input_fastas=[],
        mapping=mapping, pairs=pairs, filtered_pairs=pairs, native_admission_recheck=recheck,
        checked_records=[source, recheck])
    counts = {k: stage[k] for k in ("written_pairs", "total_pairs", "retained_pairs", "removed_mapping_pairs")}
    stage["conversion_counts"] = write("conversion_counts.json", json.dumps(counts))
    verified = {"context": stage["context"].copy(), "manifest": {"input_fastas": []}, "checked_records": []}
    if change == "context":
        stage["context"]["resolution"] = .12
    if change == "fastas":
        stage["input_fastas"] = [mapping]
    if change == "source":
        stage["source"] = pairs
    if change == "mapping":
        stage["mapping"] = pairs
    if change == "counts":
        stage["conversion_counts"] = write("conversion_counts.json", "{}")
    if change == "payload":
        write("pairs.tsv", "changed")
    write("benchmarks/results/qfo_cpm_pairs_v1/cpm_low/results.json", json.dumps(stage))
    namespaces = {"output": tmp_path / "benchmarks/results/qfo_cpm_assessment_v1/cpm_low",
                  "work": tmp_path / "qfo_benchmark/w/qcpv0",
                  "results": tmp_path / "qfo_benchmark/scoring/cpm_v1_0"}
    if change in namespaces:
        namespaces[change].mkdir(parents=True)
    if change == "symlink":
        namespaces["work"].parent.mkdir(parents=True)
        namespaces["work"].symlink_to(tmp_path / "missing")
    monkeypatch.setattr(module, "verify_sources", lambda *a: verified)
    monkeypatch.setattr(module, "environment_records", lambda *a: [])
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=0))
    monkeypatch.setattr(module.subprocess, "check_output", lambda argv, **k:
        "JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS\n22072_0|23000|COMPLETED|0:0|bizon|2\n"
        if argv[0] == "sacct" else module.CONVERTER_COMMIT)

    def command(root, converted, environment, work, results):
        assert converted == stage and environment == env
        assert work == namespaces["work"] and results == namespaces["results"]
        return ["frozen-nextflow", "six-endpoints"]

    monkeypatch.setattr(module, "command_for", command)
    if change is not None:
        with pytest.raises((ValueError, FileExistsError)):
            module.prepare(tmp_path, 0)
    else:
        prepared = module.prepare(tmp_path, 0)
        assert prepared["arm"] == "cpm_low"
        assert prepared["context"] == verified["context"]
        assert prepared["command"] == ["frozen-nextflow", "six-endpoints"]
        assert prepared["accuracy_admitted"] is False
        assert pairs in prepared["verified_records"]
        assert not namespaces["output"].exists()
