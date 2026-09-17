from copy import deepcopy

import pytest

from benchmark_tools.validate_qfo_native_assessment import AXES, swiss_families, validate_records, validate_aggregation


def records():
    rows = []
    for challenge, axes in {**AXES, "SwissTrees-X": ("TPR", "PPV")}.items():
        for metric in axes:
            rows.append({"_id": challenge + metric, "type": "assessment", "community_id": "QfO",
                         "participant_id": "fixture", "challenge_id": challenge,
                         "metrics": {"metric_id": metric, "value": 4 if metric == "NR_ORTHOLOGS" else .5, "stderr": .01}})
    return rows


@pytest.mark.parametrize("problem", [None, "missing", "duplicate", "id", "participant", "community", "unknown", "nan", "infinity", "negative", "bool", "range", "fractional_count"])
def test_native_inventory_and_numbers(problem):
    rows = records()
    if problem == "missing":
        rows.pop()
    elif problem == "duplicate":
        copy = deepcopy(rows[0]); copy["_id"] = "other"; rows.append(copy)
    elif problem == "id":
        rows[1]["_id"] = rows[0]["_id"]
    elif problem == "participant":
        rows[0]["participant_id"] = "wrong"
    elif problem == "community":
        rows[0]["community_id"] = "wrong"
    elif problem == "unknown":
        rows[0]["metrics"]["metric_id"] = "wrong"
    elif problem in {"nan", "infinity", "negative", "bool", "range"}:
        rows[0]["metrics"]["value"] = {"nan": float("nan"), "infinity": float("inf"), "negative": -.1, "bool": True, "range": 1.1}[problem]
    elif problem == "fractional_count":
        next(r for r in rows if r["metrics"]["metric_id"] == "NR_ORTHOLOGS")["metrics"]["value"] = 1.5
    if problem:
        with pytest.raises(ValueError):
            validate_records(rows, "fixture", {"X"})
    else:
        assert len(validate_records(rows, "fixture", {"X"})) == 14


@pytest.mark.parametrize("problem", [None, "axis", "value", "stderr", "participant", "extra"])
def test_native_aggregate_agreement(problem):
    native = validate_records(records(), "fixture", {"X"})
    data = {"type": "aggregation", "challenge_ids": ["GO"], "datalink": {"inline_data": {
        "visualization": {"x_axis": "NR_ORTHOLOGS", "y_axis": "avg Schlicker", "type": "2D-plot"},
        "challenge_participants": [{"participant_id": "fixture", "metric_x": 4, "metric_y": .5, "stderr_x": .01, "stderr_y": .01}]}}}
    template = deepcopy(data)
    inline = data["datalink"]["inline_data"]
    if problem == "axis":
        inline["visualization"]["y_axis"] = "F1"
    elif problem in {"value", "stderr"}:
        inline["challenge_participants"][0]["metric_y" if problem == "value" else "stderr_y"] = .4
    elif problem == "participant":
        inline["challenge_participants"][0]["participant_id"] = "wrong"
    elif problem == "extra":
        inline["challenge_participants"] *= 2
    if problem:
        with pytest.raises(ValueError):
            validate_aggregation(data, "GO", "fixture", native, template)
    else:
        assert validate_aggregation(data, "GO", "fixture", native, template)["score_semantics"] == "avg Schlicker"


@pytest.mark.parametrize("text", ["", "ReconciledTrees['X'] := RecTreeCase('Y',...)",
    "ReconciledTrees['X'] := RecTreeCase('X',...)\nReconciledTrees['X'] := RecTreeCase('X',...)"])
def test_reference_declaration_errors(text):
    with pytest.raises(ValueError):
        swiss_families(text)


def test_reference_labels_not_derived_from_prediction():
    assert swiss_families("ReconciledTrees := table():\nReconciledTrees['X'] := RecTreeCase('X',...)") == {"X"}
