from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import report_profile_evidence_pilot as readback
from benchmark_tools.metric_frontier import ENDPOINTS
from tests.unit.test_profile_evidence_readback import arm_fixture, family_cells


LABELS = ["off", "builtin_reselect", "hmmer_score_only", "hmmer_reselect"]


def cells_and_rows():
    cells = family_cells()
    cells[2] = deepcopy(cells[0])
    for i, cell in enumerate(cells):
        cell["label"] = LABELS[i]
        family = cell["details"]["family_details"]
        cell["details"].update(precision=family["precision"] / 100, recall=family["recall"] / 100)
        cell["scores"] = {ENDPOINTS[0]: family["f_score"] / 100}
        identity = "off" if i == 2 else LABELS[i]
        pin = {"path": "/" + LABELS[i], "bytes": 10, "sha256": identity}
        cell.update(products={"product.tsv": dict(pin)}, graph_stages=[{"label": "profile", "bytes": 10, "sha256": identity}],
            candidate={"candidate_partition": dict(pin), "membership_constraints": dict(pin)},
            activation={"candidate_digest": "same", "completely_measured_queries": 5, "applicability": "completely_measured"})
    control = {"f1": cells[0]["scores"][ENDPOINTS[0]], "precision": cells[0]["details"]["precision"], "recall": cells[0]["details"]["recall"]}
    rows = [{"label": label} for label in LABELS]
    rows[2]["measurement_control_matches"] = True
    for i in (1, 3):
        rows[i].update(advancement=readback.selective_gate(control, cells[i]), accuracy_interpretation="completely_measured")
    return cells, rows


def test_complete_selective_panel_verifies_measurement_identity_and_both_original_gates():
    cells, rows = cells_and_rows()
    result = readback.verify_selective_panel(cells, rows)
    assert set(result) == {"f1", "precision", "recall", "scope"}
    assert result["f1"] == 0  # Identical fixture improvements have no HMMER-specific benefit.
    assert cells[1]["advancement"]["eligible"] and cells[3]["advancement"]["eligible"]
    assert "advancement" not in cells[2]


@pytest.mark.parametrize("kind", ["missing", "order", "digest", "graph", "native", "candidate", "constraints", "scores", "families", "control_flag", "gate", "applicability"])
def test_no_drift_or_selected_rows_can_override_panel_readback(kind):
    cells, rows = cells_and_rows()
    if kind == "missing": cells.pop()
    elif kind == "order": cells.reverse()
    elif kind == "digest": cells[3]["activation"]["candidate_digest"] = "other"
    elif kind == "graph": cells[2]["graph_stages"][0]["sha256"] = "other"
    elif kind == "native": cells[2]["products"]["product.tsv"]["sha256"] = "other"
    elif kind == "candidate": cells[2]["candidate"]["candidate_partition"]["sha256"] = "other"
    elif kind == "constraints": cells[2]["candidate"]["membership_constraints"]["sha256"] = "other"
    elif kind == "scores": cells[2]["scores"][ENDPOINTS[0]] += .001
    elif kind == "families": cells[2]["details"]["family_details"]["refog_records"][0]["genes"] += 1
    elif kind == "control_flag": rows[2]["measurement_control_matches"] = False
    elif kind == "gate": rows[3]["advancement"]["eligible"] = False
    else: rows[3]["accuracy_interpretation"] = "no_complete_measurements"
    with pytest.raises(ValueError): readback.verify_selective_panel(cells, rows)


def test_incomplete_measurements_cannot_turn_point_gain_into_advancement():
    cells, rows = cells_and_rows()
    cells[3]["activation"].update(completely_measured_queries=0, applicability="no_complete_measurements")
    control = {"f1": cells[0]["scores"][ENDPOINTS[0]], "precision": cells[0]["details"]["precision"], "recall": cells[0]["details"]["recall"]}
    rows[3].update(advancement=readback.selective_gate(control, cells[3]), accuracy_interpretation="no_complete_measurements")
    readback.verify_selective_panel(cells, rows)
    assert not cells[3]["advancement"]["eligible"] and "fallback" in cells[3]["advancement"]["reason"]


def test_actual_paired_family_uncertainty_includes_backend_contrast_with_12_endpoint_tails():
    cells, _ = cells_and_rows()
    result = readback.selective_uncertainty(cells)
    assert set(result["versus_off"]["comparisons"]) == set(LABELS[1:])
    assert set(result["hmmer_versus_builtin"]["comparisons"]) == {"hmmer_reselect"}
    for key in ("versus_off", "hmmer_versus_builtin"):
        assert result[key]["replicates"] == 100000 and result[key]["seed"] == 20261018
        assert "12 planned endpoints" in result[key]["multiplicity"]
    assert "do not change point gates" in result["scope"]


@pytest.mark.parametrize("value", [True, 0, 1, 8, 9.5])
def test_multiplicity_cannot_undercover_reported_contrasts(value):
    cells, _ = cells_and_rows()
    with pytest.raises(ValueError, match="Multiplicity"):
        readback.uncertainty(cells, multiplicity_endpoints=value)


def test_selective_launch_uses_immutable_specific_snapshot_and_normalizes_plan_pin(tmp_path, monkeypatch):
    ledger = {"selective_hmmer_pilot": {"execution_plan": {"path": "plan.json", "bytes": 12, "sha256": "fixed", "input_records": 42}}}
    calls = []
    def git(command, **kw):
        calls.append(command)
        return json.dumps(ledger).encode()
    monkeypatch.setattr(readback.subprocess, "check_output", git)
    _, pin, launch = readback.launch(tmp_path, "selective_hmmer")
    assert calls == [["git", "show", readback.HMMER_LAUNCH + ":" + readback.LEDGER]]
    assert set(pin) == {"path", "bytes", "sha256"} and pin["path"] == str(tmp_path / "plan.json")
    assert launch["revision"] == readback.HMMER_LAUNCH


@pytest.mark.parametrize("mode", ["builtin_reselect", "hmmer_score_only", "hmmer_reselect"])
def test_arm_dispatch_preserves_ordered_ids_and_rechecks_native_turnover(tmp_path, monkeypatch, mode):
    activation = {"mode": mode, "candidate_digest": "same", "checked": [], "completely_measured_queries": 1}
    plan, row, native, _, folder = arm_fixture(tmp_path, monkeypatch, mode, experiment="selective_hmmer", activation=activation)
    calls = []
    def actual(graph, observed_mode, experiment, **kw):
        calls.append(kw)
        assert observed_mode == mode and experiment == "selective_hmmer"
        return activation
    monkeypatch.setattr(readback, "validate_activation", actual)
    turnover = {"candidate_pairs": 1, "added_pairs": 0, "files": [], "wall_s": 1.}
    row["native_pair_turnover"] = {**turnover, "wall_s": 2.}
    monkeypatch.setattr(readback, "native_pair_turnover", lambda *a: turnover)
    result = readback.admit_arm(tmp_path, plan, row, "selective_hmmer", ["a", "b"])
    assert calls[0]["gene_names"] == ["a", "b"] and calls[0]["folder"] == folder
    assert calls[0]["hmmer_tools"] == plan["hmmer_tools"]
    assert result["native_pair_turnover"] == turnover
    assert result["graph_stages"][0]["label"] == "fixture"
    row["native_pair_turnover"]["added_pairs"] = 1
    with pytest.raises(ValueError, match="turnover"):
        readback.admit_arm(tmp_path, plan, row, "selective_hmmer", ["a", "b"])
    with pytest.raises(ValueError, match="ordered original"):
        readback.admit_arm(tmp_path, plan, row, "selective_hmmer", {"a", "b"})
