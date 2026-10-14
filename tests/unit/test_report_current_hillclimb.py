from copy import deepcopy
import pytest

from benchmark_tools import report_current_hillclimb as reporter
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.report_current_hillclimb import upsert_cell
from benchmark_tools.run_membership_hillclimb import PRODUCTS, write_json


def test_adds_verified_endpoint_without_changing_unrelated_methods():
    untouched = {"key": "other", "label": "other", "scores": {}, "cells": []}
    methods = [deepcopy(untouched)]
    cell = {"dataset": "OrthoBench", "status": "independent_content_readback_verified",
            "scores": {"OrthoBench: OrthoBench F1": .7}}
    upsert_cell(methods, "canonical_collection:canonical", "canonical", cell)
    assert methods[0] == untouched
    assert methods[1]["scores"] == cell["scores"]
    assert methods[1]["complete_eight_endpoints"] is False


def test_replaces_stale_unscored_cell_without_inventing_scores():
    methods = [{"key": "full_search:control100", "label": "control", "scores": {},
                "cells": [{"dataset": "QfO", "status": "running", "scores": {}}]}]
    cell = {"dataset": "QfO", "status": "failed", "scores": {}, "error": "native signal"}
    upsert_cell(methods, "full_search:control100", "control", cell)
    assert methods[0]["cells"] == [cell] and methods[0]["scores"] == {}


def test_later_failed_execution_does_not_erase_verified_scores():
    good = {"dataset": "OrthoBench", "status": "independent_content_readback_verified",
            "scores": {"OrthoBench: OrthoBench F1": .7}}
    methods = [{"key": "k", "label": "k", "scores": good["scores"].copy(), "cells": [good]}]
    upsert_cell(methods, "k", "k", {"dataset": "OrthoBench", "status": "failed", "scores": {}})
    assert methods[0]["scores"] == good["scores"] and methods[0]["cells"][0] == good


def test_conflicting_verified_scores_are_not_silently_replaced():
    good = {"dataset": "OrthoBench", "status": "independent_content_readback_verified", "scores": {"a": .7}}
    methods = [{"key": "k", "label": "k", "scores": good["scores"].copy(), "cells": [good]}]
    with pytest.raises(ValueError, match="Conflicting"):
        upsert_cell(methods, "k", "k", {**good, "scores": {"a": .8}})


def native_cell():
    return {"dataset": "OrthoBench", "status": "independent_content_readback_verified",
            "scores": {"OrthoBench: OrthoBench F1": .7410607351873405},
            "products": {name: {"path": "/old/" + name, "bytes": 10, "sha256": "same"} for name in PRODUCTS}}


@pytest.mark.parametrize("name", PRODUCTS)
def test_equal_scores_cannot_hide_changed_native_predictions(name):
    good = native_cell()
    methods = [{"key": "k", "label": "k", "scores": good["scores"].copy(), "cells": [good]}]
    before = deepcopy(methods)
    changed = deepcopy(good)
    changed["products"][name]["sha256"] = "different"
    with pytest.raises(ValueError, match="despite equal scores"):
        upsert_cell(methods, "k", "k", changed)
    assert methods == before


def test_equal_scores_cannot_drop_existing_product_provenance():
    good = native_cell()
    methods = [{"key": "k", "label": "k", "scores": good["scores"].copy(), "cells": [good]}]
    changed = {name: value for name, value in good.items() if name != "products"}
    with pytest.raises(ValueError, match="missing native products"):
        upsert_cell(methods, "k", "k", changed)


def test_identical_native_products_can_move_to_another_verified_namespace():
    good = native_cell()
    methods = [{"key": "k", "label": "k", "scores": good["scores"].copy(), "cells": [good]}]
    equivalent = deepcopy(good)
    for pin in equivalent["products"].values():
        pin["path"] = pin["path"].replace("/old/", "/new/")
    upsert_cell(methods, "k", "k", equivalent)
    assert methods[0]["cells"] == [equivalent]


@pytest.mark.parametrize("status", ["failed", "not_run_failed_control", "not_run_failed_control_gate"])
def test_recovered_result_preserves_failed_and_unrun_attempts(status):
    failed = {"dataset": "OrthoBench", "status": status, "scores": {}, "producer_job_id": "original"}
    good = native_cell()
    methods = [{"key": "k", "label": "k", "scores": {}, "cells": [deepcopy(failed)]}]
    upsert_cell(methods, "k", "k", good)
    assert methods[0]["cells"] == [failed, good]
    assert methods[0]["scores"] == good["scores"]
    upsert_cell(methods, "k", "k", deepcopy(good))
    assert methods[0]["cells"] == [failed, good]


def test_repeated_failed_observation_does_not_duplicate_verified_history():
    good = native_cell()
    failed = {"dataset": "OrthoBench", "status": "failed", "scores": {}, "error": "native signal"}
    methods = [{"key": "k", "label": "k", "scores": good["scores"].copy(), "cells": [good]}]
    upsert_cell(methods, "k", "k", failed)
    upsert_cell(methods, "k", "k", deepcopy(failed))
    assert len(methods[0]["cells"]) == 2
    assert methods[0]["cells"][0] == good


def incremental_report_fixture(tmp_path, monkeypatch):
    good = native_cell()
    method = {"key": "canonical_collection:canonical", "label": "canonical",
              "scores": good["scores"].copy(), "cells": [good,
                  {"dataset": "QfO", "status": "running", "scores": {}}]}
    untouched = {"key": "canonical_collection:unregistered", "label": "historical",
                 "scores": {}, "cells": []}
    table = tmp_path / "parent.tsv"
    table.write_text("retained table\n")
    parent = tmp_path / "parent.json"
    write_json(parent, {"methods": [method, untouched], "score_table": record(table),
                        "round_screening": {}, "limitations": []})
    root = tmp_path / "canonical"
    root.mkdir()
    plan = root / "execution_plan.json"
    write_json(plan, {"inputs": [], "sources": [], "runtime": {"files": []}})
    write_json(root / "results.json", {"execution_plan": record(plan), "job_id": "live",
        "arms": [{"dataset": "QfO", "label": "canonical", "status": "independent_content_readback_verified"}]})
    new = {"dataset": "QfO", "status": "independent_content_readback_verified",
           "scores": {name: .8 for name in ENDPOINTS if name.startswith("QfO:")}}
    monkeypatch.setattr(reporter, "canonical_cell", lambda *args: deepcopy(new))
    return parent, root, method, untouched, new


def test_incremental_build_updates_existing_canonical_row(tmp_path, monkeypatch):
    parent, root, old, untouched, new = incremental_report_fixture(tmp_path, monkeypatch)
    result = reporter.build(tmp_path, parent, tmp_path / "absent_full", root, tmp_path / "out")
    canonical = next(row for row in result["methods"] if row["key"] == old["key"])
    assert canonical["scores"] == {**old["scores"], **new["scores"]}
    assert canonical["cells"][0] == old["cells"][0]
    admitted = next(row for row in canonical["cells"] if row["dataset"] == "QfO")
    assert {key: admitted[key] for key in new} == new
    assert next(row for row in result["methods"] if row["key"] == untouched["key"]) == untouched


def test_incremental_build_still_rejects_unrelated_method_changes(tmp_path, monkeypatch):
    parent, root, _, untouched, _ = incremental_report_fixture(tmp_path, monkeypatch)
    original = reporter.upsert_cell

    def corrupt_unrelated(methods, *args):
        original(methods, *args)
        next(row for row in methods if row["key"] == untouched["key"])["scores"]["invented"] = 1

    monkeypatch.setattr(reporter, "upsert_cell", corrupt_unrelated)
    with pytest.raises(ValueError, match="Unrelated method or scores changed"):
        reporter.build(tmp_path, parent, tmp_path / "absent_full", root, tmp_path / "out")
