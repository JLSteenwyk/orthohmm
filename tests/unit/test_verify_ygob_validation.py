import pytest

from benchmark_tools.benchmark_production import file_record
from benchmark_tools.verify_ygob_validation import read_metadata, require_completed_job, verify_records, verify_run


HEADER = "JobIDRaw|State|ExitCode|Elapsed\n"


def test_exact_parent_job_required():
    data = HEADER + "12|COMPLETED|0:0|01:00\n12.batch|COMPLETED|0:0|01:00\n"
    assert require_completed_job(data, 12)["JobIDRaw"] == "12"


@pytest.mark.parametrize("rows", ["", "12.batch|COMPLETED|0:0|x\n",
    "12|PENDING|0:0|x\n", "12|RUNNING|0:0|x\n", "12|FAILED|1:0|x\n",
    "12|COMPLETED|0:9|x\n", "12|COMPLETED|0:0|x\n12|COMPLETED|0:0|x\n"])
def test_bad_accounting_rejected_before_outputs_are_read(tmp_path, rows):
    with pytest.raises(ValueError, match="uniquely confirmed"):
        verify_run(tmp_path / "nonexistent", tmp_path / "also_missing", 12, HEADER + rows)


def test_artifact_bytes_and_hashes(tmp_path):
    path = tmp_path / "input.fasta"
    path.write_text(">a\nACDE\n")
    record = file_record(path, tmp_path)
    assert verify_records([record], tmp_path) == {path}
    absolute = {**record, "path": str(path)}
    assert verify_records([absolute], tmp_path) == {path}
    for invalid in ([], [record, record], [{**record, "bytes": 0}],
                    [{**record, "sha256": "0" * 64}], [{**record, "path": "../outside"}]):
        with pytest.raises(ValueError):
            verify_records(invalid, tmp_path)
    path.write_text(">a\nAAAA\n")
    with pytest.raises(ValueError, match="changed"):
        verify_records([record], tmp_path)


@pytest.mark.parametrize("text", ["job_id\t12\njob_id\t12\n", "broken\n", "\tvalue\n"])
def test_metadata_rejects_ambiguity(tmp_path, text):
    path = tmp_path / "metadata.tsv"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_metadata(path)


def test_metadata_preserves_values(tmp_path):
    path = tmp_path / "metadata.tsv"
    path.write_text("job_id\t12\nrunner_exit_code\t0\n")
    assert read_metadata(path) == {"job_id": "12", "runner_exit_code": "0"}
