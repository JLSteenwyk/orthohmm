from pathlib import Path
import csv
import json

import pytest

from benchmark_tools.assemble_wgd_application import artifact, write_table
from benchmark_tools.report_wgd_application import DIAGNOSTIC, METHODS, assemble


def admission(*rows):
    return {"native": {"checked_files": list(rows)}}


def test_artifact_accepts_identical_repeated_checks():
    row = {"path": "/native/SequenceIDs.txt", "sha256": "abc", "bytes": 10}
    assert artifact(admission(row, dict(row)), "SequenceIDs.txt") == Path(row["path"])


@pytest.mark.parametrize("change", [
    {"path": "/other/SequenceIDs.txt"}, {"sha256": "def"}, {"bytes": 11},
])
def test_artifact_rejects_conflicting_records(change):
    row = {"path": "/native/SequenceIDs.txt", "sha256": "abc", "bytes": 10}
    with pytest.raises(ValueError, match="Missing or ambiguous"):
        artifact(admission(row, {**row, **change}), "SequenceIDs.txt")


def test_artifact_rejects_missing_exact_basename():
    with pytest.raises(ValueError, match="Missing or ambiguous"):
        artifact(admission({"path": "/native/old_SequenceIDs.txt"}), "SequenceIDs.txt")


def test_table_keeps_exclusion_and_all_method_fields(tmp_path):
    owners = {"a": "Scerevisiae", "b": "Scerevisiae", "c": "Smikatae"}
    cohort = [
        {"orf_pair": ["a", "b"], "experimental_class": "High", "split_eligible": True,
         "reference_eligible": True, "reference_pillar": "p", "available_pillar_members": list(owners),
         "available_members_by_species": {"Scerevisiae": 2, "Smikatae": 1}},
        {"orf_pair": ["a", "z"], "experimental_class": "Low", "split_eligible": False,
         "reference_eligible": False, "input_reasons": ["absent"]},
    ]
    admitted = {m: {"status": "admitted", "groups": {"x": ["a", "c"], "y": ["b"]}}
                for m in (*METHODS, DIAGNOSTIC)}
    report = assemble({"cohort_pairs": cohort, "prespecified_examples": []}, {"p": list(owners)}, owners, admitted)
    path = tmp_path / "pairs.tsv"
    write_table(report, path)
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 2
    assert len(rows[0]) == 76
    assert json.loads(rows[1]["input_reasons"]) == ["absent"]
    for method in (*METHODS, DIAGNOSTIC):
        assert rows[0][method + ".separation_rate"] == "1"
        assert rows[0][method + ".supported_separation_rate"] == "0"
        assert rows[1][method + ".assignment_state"] == "input_excluded"
        assert rows[1][method + ".separation_rate"] == ""
    with pytest.raises(FileExistsError):
        write_table(report, path)
