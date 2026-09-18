import json

import pytest

from benchmark_tools import audit_orthomcl_search_table as module


def fixture(tmp_path):
    fasta, blast, log = (tmp_path / name for name in ("all.fa", "all.blast", "blast.log"))
    fasta.write_text(">A\nACDEFGHIKL\n>B\nACDEFGHIKL\n>C\nXXXX\n>D\nAAAA\n")
    blast.write_text("A\tA\t100.00\t10\t0\t0\t1\t10\t1\t10\t1e-9\t40\n"
                     "A\tC\t100.00\t4\t0\t0\t1\t4\t1\t4\te-10\t30\n")
    log.write_text("[blastall] WARNING: [000.000] C: SetUpBlastSearch failed.\n")
    return blast, fasta, log


def test_query_failure_is_distinct_from_incoming_hits_and_no_hits(tmp_path):
    content = module.audit_table(*fixture(tmp_path))
    assert content["hsp_rows"] == content["distinct_directed_pairs"] == 2
    assert content["queries_with_hits"] == 1
    assert content["subjects_with_hits"] == 2
    assert content["failed_queries_with_subject_hits"] == 1
    assert content["failed_queries_with_query_hits"] == 0
    assert content["query_ids_without_hits"] == ["B", "C", "D"]
    assert content["queries_without_hits_and_without_logged_failure"] == 2


@pytest.mark.parametrize("column,value", [(0, "unknown"), (1, "unknown"), (2, "nan"),
    (2, "101"), (2, "50"), (3, "0"), (3, "10.5"), (4, "-1"), (4, "11"),
    (5, "1"), (6, "0"), (7, "11"), (8, "10"), (9, "1"),
    (10, "nan"), (10, "-1"), (10, "inf"), (11, "-1"), (11, "nan")])
def test_bad_fields_rejected(tmp_path, column, value):
    paths = fixture(tmp_path)
    fields = paths[0].read_text().splitlines()[0].split("\t")
    fields[column] = value
    paths[0].write_text("\t".join(fields) + "\n")
    with pytest.raises(ValueError):
        module.audit_table(*paths)


@pytest.mark.parametrize("text", ["", "\n", "A\tB\n"])
def test_empty_or_malformed_table(tmp_path, text):
    paths = fixture(tmp_path)
    paths[0].write_text(text)
    with pytest.raises(ValueError):
        module.audit_table(*paths)


def test_contiguous_hsps_and_gaps(tmp_path):
    paths = fixture(tmp_path)
    paths[0].write_text("A\tB\t80.00\t10\t1\t1\t1\t9\t1\t10\t1e-8\t40\n" * 2)
    content = module.audit_table(*paths)
    assert content["hsp_rows"] == 2 and content["distinct_directed_pairs"] == 1


def test_rounded_identity_and_above_cutoff_hsp_are_preserved(tmp_path):
    paths = fixture(tmp_path)
    paths[0].write_text("A\tB\t66.67\t3\t1\t0\t1\t3\t1\t3\t2e-5\t4\n")
    content = module.audit_table(*paths)
    assert content["hsp_rows"] == content["hsp_rows_above_1e_minus_5"] == 1


@pytest.mark.parametrize("kind", ["query", "subject"])
def test_noncontiguous_blocks_rejected(tmp_path, kind):
    paths = fixture(tmp_path)
    line = paths[0].read_text().splitlines()[0] + "\n"
    other = line.replace("A\tA", "B\tA" if kind == "query" else "A\tB")
    paths[0].write_text(line + other + line)
    with pytest.raises(ValueError, match="Noncontiguous " + kind):
        module.audit_table(*paths)


def test_unknown_diagnostic_and_empty_sequence_rejected(tmp_path):
    paths = fixture(tmp_path)
    paths[2].write_text("[blastall] WARNING: [0] missing: SetUpBlastSearch failed.\n")
    with pytest.raises(ValueError, match="absent"):
        module.audit_table(*paths)
    paths[1].write_text(">empty\n")
    with pytest.raises(ValueError, match="nonempty"):
        module.audit_table(*paths)


def test_report_hashes_and_no_overwrite(tmp_path):
    paths = fixture(tmp_path)
    output = tmp_path / "report.json"
    report = module.audit(*paths, output)
    assert json.loads(output.read_text()) == report
    assert report["search_admitted"] is False and report["accuracy_admitted"] is False
    for item in report["checked_records"]:
        module.check(item)
    with pytest.raises(FileExistsError):
        module.audit(*paths, output)


def test_changed_input_during_audit_rejected(tmp_path, monkeypatch):
    paths = fixture(tmp_path)
    original = module.audit_table
    def changed(*args):
        result = original(*args)
        paths[2].write_text("changed")
        return result
    monkeypatch.setattr(module, "audit_table", changed)
    output = tmp_path / "report.json"
    with pytest.raises(ValueError):
        module.audit(*paths, output)
    assert not output.exists()
