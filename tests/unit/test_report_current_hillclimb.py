from copy import deepcopy
import pytest

from benchmark_tools.report_current_hillclimb import upsert_cell
from benchmark_tools.run_membership_hillclimb import PRODUCTS


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
