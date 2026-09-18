import pytest

from benchmark_tools import audit_orthomcl_bpo_content as module
from benchmark_tools.probe_orthomcl_bpo_parity import fixture
from benchmark_tools.convert_orthomcl_blast import convert_blast


def paths(tmp_path):
    fasta, blast = fixture(tmp_path)
    bpo = tmp_path / "all.bpo"
    convert_blast(blast, fasta, bpo, progress_every=0)
    return blast, fasta, bpo


def test_complete_content_with_multi_hsp_gap_and_cutoff(tmp_path):
    assert module.validate(*paths(tmp_path)) == {
        "source_hsp_rows": 8, "source_pair_blocks": 7, "cutoff_excluded_pair_blocks": 1,
        "bpo_pair_records": 6, "input_proteins": 4}


def test_first_hsp_controls_pair_cutoff_not_minimum_hsp(tmp_path):
    args = paths(tmp_path)
    args[0].write_text(args[0].read_text().replace("e-20\t100", "2e-5\t100"))
    convert_blast(args[0], args[1], args[2], progress_every=0)
    counts = module.validate(*args)
    assert counts["cutoff_excluded_pair_blocks"] == 2
    assert counts["bpo_pair_records"] == 5


def test_second_hsp_span_is_checked(tmp_path):
    args = paths(tmp_path)
    args[0].write_text(args[0].read_text().replace("11\t18\t3\t8", "12\t19\t3\t8"))
    with pytest.raises(ValueError, match="fields differ"):
        module.validate(*args)


@pytest.mark.parametrize("column,value", [(0, "2"), (1, "B"), (2, "119"), (3, "D"),
    (4, "119"), (5, "1e-4"), (6, "99"), (7, "1:1-119:1-120.")])
def test_each_bpo_field_is_checked(tmp_path, column, value):
    args = paths(tmp_path)
    lines = args[2].read_text().splitlines()
    fields = lines[0].split(";")
    fields[column] = value
    lines[0] = ";".join(fields)
    args[2].write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="fields differ"):
        module.validate(*args)


@pytest.mark.parametrize("kind", ["extra", "missing", "empty", "swap"])
def test_bpo_cardinality_and_order(tmp_path, kind):
    args = paths(tmp_path)
    lines = args[2].read_text().splitlines(keepends=True)
    if kind == "extra":
        lines.append(lines[-1])
    elif kind == "missing":
        lines.pop()
    elif kind == "empty":
        lines = []
    else:
        lines[:2] = reversed(lines[:2])
    args[2].write_text("".join(lines))
    with pytest.raises(ValueError):
        module.validate(*args)


@pytest.mark.parametrize("kind", ["query_repeat", "subject_repeat", "unknown", "nan", "coordinates", "gap", "columns"])
def test_source_corruption(tmp_path, kind):
    args = paths(tmp_path)
    lines = args[0].read_text().splitlines(keepends=True)
    if kind == "query_repeat":
        lines.append(lines[0])
    elif kind == "subject_repeat":
        lines.insert(4, lines[0])
    else:
        fields = lines[0].strip().split("\t")
        if kind == "unknown":
            fields[0] = "unknown"
        elif kind == "nan":
            fields[10] = "NaN"
        elif kind == "coordinates":
            fields[7] = "121"
        elif kind == "gap":
            fields[5] = "1"
        else:
            fields.pop()
        lines[0] = "\t".join(fields) + "\n"
    args[0].write_text("".join(lines))
    with pytest.raises(ValueError):
        module.validate(*args)


def test_report_records_and_no_overwrite(tmp_path):
    args = paths(tmp_path)
    output = tmp_path / "report.json"
    result = module.audit(*args, output)
    assert result["status"] == "bpo_content_matches_source_hsps"
    assert result["accuracy_admitted"] is False
    for item in result["checked_records"]:
        module.check(item)
    with pytest.raises(FileExistsError):
        module.audit(*args, output)


def test_mutation_during_validation_rejected(tmp_path, monkeypatch):
    args = paths(tmp_path)
    original = module.validate
    def changed(*values):
        content = original(*values)
        args[2].write_text("changed")
        return content
    monkeypatch.setattr(module, "validate", changed)
    output = tmp_path / "report.json"
    with pytest.raises(ValueError):
        module.audit(*args, output)
    assert not output.exists()
