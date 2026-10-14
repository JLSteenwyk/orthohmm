from copy import deepcopy
import pytest

from benchmark_tools.report_current_hillclimb import upsert_cell


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
