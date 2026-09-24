import gzip

import pytest

from benchmark_tools.audit_orthomcl_reference_impact import audit_vgnc, parse_native_report, failed_query_targets


REPORT = """Darwin banner
QFOAUDIT\tFAMILY\tSwissTrees\tS1\t10
QFOAUDIT\tMEMBER\tSwissTrees\tS1\t5
QFOAUDIT\tRELATION\tSwissTrees\tS1\t5\t6\tS
QFOAUDIT\tFAMILY\tTreeFam-A\tT1\t4
QFOAUDIT\tANNOTATION\t5\t0\t2
QFOAUDIT\tCOMPLETE\t1
"""


def test_parse_native_report():
    families, annotations = parse_native_report(REPORT, {5})
    assert len(families) == 2
    assert families[0]["incident_relations"] == [{"pair": [5, 6], "event": "S"}]
    assert annotations == {5: {"ec_terms": 0, "experimental_go_terms": 2}}


@pytest.mark.parametrize("text", [
    REPORT.replace("QFOAUDIT\tCOMPLETE\t1\n", ""),
    REPORT.replace("ANNOTATION\t5", "ANNOTATION\t7"),
    REPORT.replace("\t5\t6\tS", "\t6\t5\tS"),
    REPORT.replace("\t5\t6\tS", "\t6\t7\tS"),
    REPORT.replace("\t5\t0\t2", "\t5\t-1\t2"),
    REPORT + "QFOAUDIT\tCOMPLETE\t1\n",
    REPORT + "Error, interrupted tree read\n",
])
def test_reject_invalid_native_report(text):
    with pytest.raises(ValueError):
        parse_native_report(text, {5})


def test_vgnc_incident_pairs_are_unique(tmp_path):
    path = tmp_path / "vgnc.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("1\t2\tA\n2\t1\tA\n1\t3\tA\n4\t5\tB\n")
    result = audit_vgnc(path, {1, 2})
    assert result["reference_pairs"] == 3
    assert len(result["incident_pairs"]) == 2
    assert result["incident_pair_fraction"] == pytest.approx(2 / 3)
    with gzip.open(path, "wt") as handle:
        handle.write("1\t2\tA\n2\t1\tB\n")
    with pytest.raises(ValueError, match="Conflicting"):
        audit_vgnc(path, {1})


def target_fixture():
    return {"failed_queries": 2, "records": [
        dict(gene="sp|A|a", accession="A", query_failed=True, final_group_line=None, final_group_size=0),
        dict(gene="sp|B|b", accession="B", query_failed=True, final_group_line=None, final_group_size=0)]}


def test_failed_query_targets_default_and_grouped():
    data = target_fixture()
    assert failed_query_targets(data, {"A": 1, "B": 2})[1] == {1, 2}
    data["records"][1].update(final_group_line=3, final_group_size=5)
    with pytest.raises(ValueError, match="opt in"):
        failed_query_targets(data, {"A": 1, "B": 2})
    failed, targets = failed_query_targets(data, {"A": 1, "B": 2}, allow_grouped=True)
    assert targets == {1, 2}
    assert [r["final_group_line"] for r in failed] == [None, 3]


@pytest.mark.parametrize("problem", ["count", "bool_count", "flag", "duplicate_gene", "duplicate_accession",
    "mapping_collision", "bool_mapping", "empty", "ungrouped_size", "grouped_size", "line"])
def test_failed_query_target_rejection(problem):
    data, mapping = target_fixture(), {"A": 1, "B": 2}
    if problem == "count":
        data["failed_queries"] = 1
    elif problem == "bool_count":
        data["failed_queries"] = True
    elif problem == "flag":
        data["records"][0]["query_failed"] = 1
    elif problem == "duplicate_gene":
        data["records"][1]["gene"] = data["records"][0]["gene"]
    elif problem == "duplicate_accession":
        data["records"][1]["accession"] = "A"
    elif problem == "mapping_collision":
        mapping["B"] = 1
    elif problem == "bool_mapping":
        mapping["A"] = True
    elif problem == "empty":
        data = dict(failed_queries=0, records=[])
    elif problem == "ungrouped_size":
        data["records"][0]["final_group_size"] = 2
    elif problem == "grouped_size":
        data["records"][0]["final_group_line"] = 1
    elif problem == "line":
        data["records"][0].update(final_group_line=True, final_group_size=2)
    with pytest.raises(ValueError):
        failed_query_targets(data, mapping, allow_grouped=True)
