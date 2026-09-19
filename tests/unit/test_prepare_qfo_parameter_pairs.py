import gzip
import json
from pathlib import Path

import pytest

from benchmark_tools import prepare_qfo_parameter_pairs as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture():
    pairs = {"path": "native.tsv", "sha256": "abc", "bytes": 10}
    return {"status": "corrected_qfo_parameter_native_pairs_verified", "variant": "norm_low", "index": 0,
        "cell": {"label": "candidate_norm_low"}, "source": {"path": "validator"},
        "accuracy_evaluated": False, "scoring_admitted": False, "publication_ready": False,
        "native_pair_count": 2, "native_pairs": pairs, "checked_records": [pairs]}


def test_valid_native_identity():
    module.validate_native(fixture(), 0, {"path": "validator"})


@pytest.mark.parametrize("key,value", [("status", "prepared"), ("variant", "margin_low"), ("index", 1),
    ("cell", {"label": "p1_c1_r1"}), ("source", {}), ("accuracy_evaluated", True),
    ("scoring_admitted", True), ("publication_ready", True), ("native_pair_count", True),
    ("native_pair_count", 0), ("checked_records", [])])
def test_wrong_admission_rejected(key, value):
    native = fixture()
    native[key] = value
    with pytest.raises(ValueError):
        module.validate_native(native, 0, {"path": "validator"})


@pytest.mark.parametrize("index", [-1, 4, True, "0", 0.0])
def test_invalid_index_rejected(index):
    with pytest.raises(ValueError):
        module.validate_native({}, index, {})
    with pytest.raises(ValueError):
        module.completed_admission("", index)


@pytest.mark.parametrize("state,cpus,exit_code", [("RUNNING", "2", "0:0"), ("FAILED", "2", "1:0"),
                                               ("COMPLETED", "32", "0:0"), ("COMPLETED", "2", "1:0")])
def test_scheduler_gate(state, cpus, exit_code):
    accounting = ("JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS\n"
                  f"21935_0|21936|{state}|{exit_code}|bizon|{cpus}\n")
    with pytest.raises(ValueError):
        module.completed_admission(accounting, 0)


def test_scheduler_selects_array_not_raw_id():
    accounting = "JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS\n21935_0|21936|COMPLETED|0:0|bizon|2\n"
    assert module.completed_admission(accounting, 0)["JobIDRaw"] == "21936"


@pytest.mark.parametrize("missing", [False, True])
def test_actual_native_conversion_preserves_pairs_and_records_mapping_loss(tmp_path, missing):
    native_file = tmp_path / "native.tsv"
    native_file.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\nA\ts1\tB\ts2\nA\ts1\tC\ts3\n")
    native = {"native_pairs": {"path": str(native_file)}, "native_pair_count": 2}
    mapping = tmp_path / "mapping.json.gz"
    with gzip.open(mapping, "wt") as stream:
        json.dump({"mapping": dict.fromkeys("AB" if missing else "ABC", 1)}, stream)
    owners = {"A": "s1", "B": "s2", "C": "s3"}
    if missing:
        with pytest.raises(ValueError, match="mapping loss"):
            module.convert(native, owners, mapping, tmp_path)
    else:
        counts = module.convert(native, owners, mapping, tmp_path)
        assert counts["written_pairs"] == counts["retained_pairs"] == 2
    counts = json.loads((tmp_path / "conversion_counts.json").read_text())
    assert counts["removed_mapping_pairs"] == int(missing)
    assert (tmp_path / "pairs.partial.tsv").read_text() == "A\tB\nA\tC\n"
    assert (tmp_path / "pairs.qfo.partial.tsv").read_text() == ("A\tB\n" if missing else "A\tB\nA\tC\n")
    assert not (tmp_path / "pairs.tsv").exists()


def test_allocation_rejected_before_reading_admission(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: pytest.fail("Read scheduler"))
    with pytest.raises(ValueError, match="scheduled"):
        module.prepare(tmp_path, 0)


@pytest.mark.parametrize("mismatch", [False, True])
def test_prepare_with_real_conversion_and_mocked_frozen_recheck(tmp_path, monkeypatch, mismatch):
    executor = tmp_path / "benchmarks/work/publication_qfo_parameter_native_admission_v1"
    checker = executor / "benchmark_tools/admit_qfo_parameter_phylogeny.py"
    checker.parent.mkdir(parents=True)
    checker.write_text("# frozen fixture\n")
    pairs = tmp_path / "native.tsv"
    pairs.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\nA\ts1\tB\ts2\nA\ts1\tC\ts3\n")
    metadata = tmp_path / "metadata.json"
    metadata.write_text("{}")
    native = fixture()
    native.update(source=record(checker), native_pairs=record(pairs), checked_records=[record(pairs)],
                  helpers=[], native_group_integrity={"native_manifest": record(metadata)})
    admission = tmp_path / "benchmarks/work/qfo_parameter_native_admission_21935_0.json"
    admission.write_text(json.dumps(native))
    mapping = tmp_path / "mapping.json.gz"
    with gzip.open(mapping, "wt") as stream:
        json.dump({"mapping": dict.fromkeys("ABC", 1)}, stream)
    env_path = tmp_path / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    env_path.parent.mkdir(parents=True)
    env_path.write_text(json.dumps({"reference_files": [record(mapping)]}))
    monkeypatch.setattr(module, "read_frozen", lambda path, sha: json.loads(path.read_text()))
    monkeypatch.setattr(module, "verify_sources", lambda *a: (
        {"partition": {"path": "candidate"}}, {"input_fastas": []}, None, None, None, None, [], None))
    monkeypatch.setattr(module, "gene_ownership", lambda *a: ({"A": "s1", "B": "s2", "C": "s3"}, {}))
    for key, value in {"SLURM_JOB_ID": "23100", "SLURM_CPUS_PER_TASK": "2",
                       "SLURM_JOB_NODELIST": "bizon", "SLURM_ARRAY_TASK_ID": "0"}.items():
        monkeypatch.setenv(key, value)
    def check_output(command, **kwargs):
        if command[0] == "sacct":
            return "JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS\n21935_0|21936|COMPLETED|0:0|bizon|2\n"
        return module.ADMITTER_COMMIT + "\n"
    monkeypatch.setattr(module.subprocess, "check_output", check_output)
    def run(command, **kwargs):
        if command[0] != "git":
            assert command[1] == str(checker)
            assert command[command.index("--index") + 1] == "0"
            Path(command[command.index("--output") + 1]).write_text(json.dumps({} if mismatch else native))
    monkeypatch.setattr(module.subprocess, "run", run)
    output = tmp_path / "benchmarks/results/qfo_parameter_pairs_v1/norm_low"
    if mismatch:
        with pytest.raises(ValueError, match="Fresh native"):
            module.prepare(tmp_path, 0)
        assert json.loads((output / "results.json").read_text())["status"] == "failed"
        assert not (output / "pairs.tsv").exists()
    else:
        report = module.prepare(tmp_path, 0)
        assert report["status"] == "corrected_parameter_native_pairs_prepared_unscored"
        assert report["native_admission_recheck"] in report["checked_records"]
        assert (output / "pairs.qfo.tsv").read_text() == "A\tB\nA\tC\n"
        assert report["accuracy_evaluated"] is False
    with pytest.raises(FileExistsError):
        module.prepare(tmp_path, 0)
