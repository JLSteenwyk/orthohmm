from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools.audit_publication_figures import validate_ob_complete


def manifest():
    root = Path(__file__).resolve().parents[2]
    return json.loads((root / "benchmark_tools/results/ob_complete_strata_figure_20260928/manifest.json").read_text())


def test_complete_panel():
    validate_ob_complete(manifest())


@pytest.mark.parametrize("key,value", [("displayed_finite_values", 287), ("missing_values", 0), ("new_intervals", True)])
def test_reject_partial_or_inferential_panel(key, value):
    data = deepcopy(manifest())
    data[key] = value
    with pytest.raises(ValueError):
        validate_ob_complete(data)


def test_reject_unpinned_source():
    data = manifest()
    data["source"]["sha256"] = "0" * 64
    with pytest.raises(ValueError):
        validate_ob_complete(data)
