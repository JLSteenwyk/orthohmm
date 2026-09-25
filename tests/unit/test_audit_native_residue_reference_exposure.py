import json
from pathlib import Path

import pytest

from benchmark_tools import audit_native_residue_reference_exposure as module


def fixture():
    review = json.loads((Path(module.__file__).parent /
        "results/qfo_recovery_native_residue_deletions_20260925.json").read_text())
    mapping = {t["id"].split("|")[1]: i for i, t in enumerate(review["transformations"], 1)}
    owners = {t["id"]: "taxon" for t in review["transformations"]}
    return review, mapping, owners


@pytest.mark.parametrize("problem", ["duplicate", "number", "identity", "unknown", "residue"])
def test_bad_selection_rejected(problem):
    review, mapping, owners = fixture()
    first = review["transformations"][0]
    accession = first["id"].split("|")[1]
    if problem == "duplicate":
        mapping[accession] = 2
    elif problem == "number":
        mapping[accession] = True
    elif problem == "identity":
        first["id"] = "malformed"
    elif problem == "unknown":
        owners.clear()
    else:
        first["deleted_residue"] = "X"
    with pytest.raises(ValueError):
        module.select_targets(review, mapping, owners)


def test_reference_summary_does_not_infer_search_failure_or_score_effect():
    rows = module.select_targets(*fixture())
    lines = ["QFOAUDIT\tFAMILY\tSwissTrees\tSwiss1\t6",
             "QFOAUDIT\tFAMILY\tTreeFam-A\tTF1\t5",
             "QFOAUDIT\tMEMBER\tSwissTrees\tSwiss1\t1",
             "QFOAUDIT\tMEMBER\tTreeFam-A\tTF1\t2",
             "QFOAUDIT\tRELATION\tSwissTrees\tSwiss1\t1\t8\tS"]
    lines += [f"QFOAUDIT\tANNOTATION\t{i}\t1\t0" for i in range(1, 8)]
    lines += ["QFOAUDIT\tCOMPLETE\t7"]
    fas = {"taxon": {"feature": {rows[0]["accession"]: {"pfam": {"domain": {}}, "length": 10}}}}
    result = module.summarize(rows, "\n".join(lines), {"incident_pairs": []}, fas)
    assert result["selected_proteins"] == 7
    assert result["tree_summary"]["SwissTrees"]["eligible_cases_with_selected_members"] == 1
    assert result["tree_summary"]["TreeFam-A"]["eligible_cases"] == 0
    assert result["records"][0]["fas_feature_types_by_tool"] == {"pfam": 1}
    assert result["records"][1]["fas_annotation_entry_present"] is False
    assert "failed_members" not in json.dumps(result)
    assert "query_failed" not in json.dumps(result)
    assert result["accuracy_admitted"] is False and result["exact_sequence_parity"] is False
    with pytest.raises(ValueError, match="Incomplete"):
        module.summarize(rows, "\n".join(lines[:-1]), {}, fas)


def test_existing_output_rejected(tmp_path):
    with pytest.raises(FileExistsError):
        module.run(tmp_path, tmp_path)
