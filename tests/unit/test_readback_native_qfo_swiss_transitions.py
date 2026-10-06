"""Actual-artifact SQL readback and deliberately corrupted diagnostic reports."""

import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import readback_native_qfo_swiss_transitions as sql

RESULTS = Path(sql.__file__).parent / "results"
REPORT = RESULTS / "native_qfo_swiss_pair_transitions_20261006_v2.json"
REPORT_SHA = "fe003f2cbc4285ea56cd80a92b71703b5c9e135ad1c8f244504dcfeab0911b46"


def test_actual_artifact_independent_sql_readback():
    result = sql.verify(REPORT, REPORT_SHA)
    assert (result["families_checked"], result["raw_rows_checked"], result["paired_relations_checked"],
            result["changed_relations_checked"], result["transition_cells_checked"]) == (18, 21530, 10765, 2023, 288)
    retained = json.loads((RESULTS / "native_qfo_swiss_pair_transition_sql_readback_20261006.json").read_text())
    assert result == retained
    assert result["uncertainty_admitted"] is result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ("family_marginal", "family_table", "pooled_table", "family_summary", "pooled_summary",
    "subset", "macro", "native_endpoint", "timing", "scope", "source", "missing_ledger_row", "duplicate_ledger_row",
    "unrecorded_raw", "raw_hash", "family_inventory", "reference_count", "changed_count", "cell_order"))
def test_corrupt_actual_report_refuses_before_success(tmp_path, fault):
    report = copy.deepcopy(json.loads(REPORT.read_text()))
    c = report["comparison"]
    if fault == "family_marginal":
        c["families"][0]["before_counts"]["TP"] += 1
    elif fault == "family_table":
        c["families"][0]["transitions"]["TP->FN"] += 1
    elif fault == "pooled_table":
        c["transitions"]["TP->FN"] += 1
    elif fault == "family_summary":
        c["families"][0]["removed_true_positives"] += 1
    elif fault == "pooled_summary":
        c["removed_true_positives"] += 1
    elif fault == "subset":
        c["after_predictions_subset_on_reference"] = False
    elif fault == "macro":
        report["cells"][0]["macro_statistics"]["F1"] += .01
    elif fault == "native_endpoint":
        report["cells"][0]["native_endpoint_f1"] += .01
    elif fault == "timing":
        report["cells"][1]["timing_eligible"] = True
    elif fault == "scope":
        report["uncertainty_admitted"] = True
    elif fault == "source":
        report["source"] = sql.record(__file__)
    elif fault in ("missing_ledger_row", "duplicate_ledger_row"):
        lines = Path(report["changed_relations_ledger"]["path"]).read_text().splitlines()
        lines = lines[:-1] if fault == "missing_ledger_row" else lines + [lines[-1]]
        ledger = tmp_path / "changes.tsv"
        ledger.write_text("\n".join(lines) + "\n")
        report["changed_relations_ledger"] = sql.record(ledger)
    elif fault == "unrecorded_raw":
        report["checked_inputs"] = [ref for ref in report["checked_inputs"] if ref != report["cells"][0]["raw"]]
    elif fault == "raw_hash":
        report["checked_inputs"][0]["sha256"] = "0" * 64
    elif fault == "family_inventory":
        c["families"].pop()
    elif fault == "reference_count":
        c["reference_relations"] += 1
    elif fault == "changed_count":
        c["changed_relations"] += 1
    elif fault == "cell_order":
        report["cells"].reverse()
    path = tmp_path / "report.json"
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError):
        sql.verify(path, sql.record(path)["sha256"])


def test_changed_report_digest_refuses():
    with pytest.raises(ValueError, match="Changed diagnostic report"):
        sql.verify(REPORT, "0" * 64)
