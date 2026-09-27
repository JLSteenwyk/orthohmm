import copy

import pytest

from benchmark_tools.admit_qfo_raw_tree_reuse import (
    ARCHIVED_ORCHESTRATION, audited_records, collect_records, record, verify_completion,
)


def test_nested_records_deduplicated_and_normalized():
    r = dict(path="/x", bytes=3, sha256="abc")
    assert collect_records([r, dict(nested=[dict(r, label="extra")])]) == [r]


def test_conflicting_records_rejected():
    r = dict(path="/x", bytes=3, sha256="abc")
    with pytest.raises(ValueError):
        collect_records([r, dict(r, sha256="changed")])


def binding():
    plan, result_record, admission = [dict(path="/" + n, bytes=1, sha256=n) for n in ("plan", "result", "admission")]
    execution = dict(status="readback_complete", job_id="22332", plan=plan, result=result_record)
    result = dict(status="fresh_qfo_phylogeny_scientific_readback_complete", job_id=22329,
                  admission=admission, accuracy_evaluated=False, historical_scores_replaced=False,
                  partition=dict(label_invariant_equal=False))
    return execution, result, plan, result_record, admission


def test_valid_differences_do_not_block_reuse():
    verify_completion(*binding())


@pytest.mark.parametrize("index,key,value", [(0,"status","partial"), (0,"job_id","22329"), (0,"job_id","22330"),
    (0,"plan",{}), (0,"result",{}), (1,"status","partial"), (1,"job_id",22328),
    (1,"admission",{}), (1,"accuracy_evaluated",True), (1,"historical_scores_replaced",True)])
def test_completion_mutations_rejected(index, key, value):
    args = copy.deepcopy(binding())
    args[index][key] = value
    with pytest.raises(ValueError):
        verify_completion(*args)


@pytest.mark.parametrize("mutation", [None, "archive", "scientific_source", "missing_identity"])
def test_audited_archive_is_narrow_and_content_bound(tmp_path, mutation):
    repo, directory = tmp_path / "repo", tmp_path / "run"
    (repo / "benchmark_tools").mkdir(parents=True)
    archive = directory / "readback_v2_source_archive"
    archive.mkdir(parents=True)
    records = []
    for name in ARCHIVED_ORCHESTRATION:
        original = repo / "benchmark_tools" / name
        original.write_text("audited")
        records.append(record(original))
        (archive / name).write_text("audited")
        original.write_text("updated orchestration")
    science = repo / "benchmark_tools/scientific_reader.py"
    science.write_text("frozen")
    records.append(record(science))
    if mutation == "archive":
        (archive / ARCHIVED_ORCHESTRATION[0]).write_text("changed")
    elif mutation == "scientific_source":
        science.write_text("changed")
    elif mutation == "missing_identity":
        records.pop(0)
    if mutation:
        with pytest.raises(ValueError):
            audited_records(repo, directory, records)
    else:
        verified, mapping = audited_records(repo, directory, records)
        assert len(mapping) == 3
        assert verified[-1] == record(science)
        assert all(item["original"]["sha256"] == item["archived"]["sha256"] for item in mapping)
