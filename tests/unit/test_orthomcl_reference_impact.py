import gzip

import pytest

from benchmark_tools.audit_orthomcl_reference_impact import audit_vgnc, parse_native_report


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
