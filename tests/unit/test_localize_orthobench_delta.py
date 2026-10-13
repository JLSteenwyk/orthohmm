from copy import deepcopy
import json

import pytest

from benchmark_tools import localize_orthobench_delta as helper
from benchmark_tools.localize_orthobench_delta import localize, reconstruct
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.score_orthobench_partition import score_partition


def fixture():
    references = {"family1": {"a", "b", "c"}, "family2": {"d", "e"}}
    predictions = [[{"a", "b", "c"}, {"d", "e"}, {"x"}], [{"a", "b", "c", "x"}, {"d", "e"}]]
    cells = []
    for groups in predictions:
        details = score_partition(groups, references, {})
        cells.append({"dataset": "OrthoBench", "status": "independent_content_readback_verified",
            "details": {"family_details": details}, "scores": {"OrthoBench: OrthoBench F1": details["f_score"] / 100}})
    return {"cells": cells}


def test_localization_reproduces_weighted_endpoints_without_additive_f1_claim():
    result = localize(fixture())
    assert result["families"] == 2 and result["changed_families"] == 1
    assert result["weighted_delta"] == {"true_positive": 0, "false_positive": 1.5, "false_negative": 0}
    assert result["control"]["F1"] == 1 and result["variant"]["F1"] < 1
    assert result["changes"][0]["refog"] == "family1" and result["defaults_promoted"] is False


@pytest.mark.parametrize("field", ["precision", "recall", "f_score", "weighted_true_positive",
                                   "weighted_false_positive", "weighted_false_negative"])
def test_reconstruction_rejects_inconsistent_readback(field):
    details = fixture()["cells"][0]["details"]["family_details"]
    details[field] += .1
    with pytest.raises(ValueError):
        reconstruct(details)


def test_different_reference_universe_is_not_a_comparable_experiment():
    summary = fixture()
    summary["cells"][1]["details"]["family_details"]["refog_records"][0]["refog"] = "different"
    with pytest.raises(ValueError, match="Reference families"):
        localize(summary)


def test_unverified_or_nonorthobench_cell_is_not_admitted():
    for key, value in (("status", "scored_development_result"), ("dataset", "QfO")):
        summary = fixture()
        summary["cells"][1][key] = value
        with pytest.raises(ValueError, match="independently verified"):
            localize(summary)


def test_duplicate_family_and_unreproduced_claim_are_rejected():
    summary = fixture()
    details = summary["cells"][0]["details"]["family_details"]
    details["refog_records"].append(deepcopy(details["refog_records"][0]))
    with pytest.raises(ValueError, match="Duplicate family"):
        localize(summary)
    summary = fixture()
    summary["cells"][0]["scores"]["OrthoBench: OrthoBench F1"] = .5
    with pytest.raises(ValueError, match="Claimed cell F1"):
        localize(summary)


def test_cli_resolves_compact_summary_to_pinned_full_readback(tmp_path, monkeypatch):
    summary = fixture()
    full = deepcopy(summary)
    for cell in full["cells"]:
        cell["checked"] = []
    path = tmp_path / "full.json"
    path.write_text(json.dumps(full))
    summary["report"] = record(path)
    compact = tmp_path / "summary.json"
    compact.write_text(json.dumps(summary))
    output = tmp_path / "localization.json"
    monkeypatch.setattr(helper.sys, "argv", ["localizer", "--summary", str(compact), "--output", str(output)])
    helper.main()
    result = json.loads(output.read_text())
    assert result["full_readback"] == record(path) and result["input_pins_reverified"] is True
    assert result["changed_families"] == 1
