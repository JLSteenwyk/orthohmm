from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import reproduce_orthobench_factorial as module


def fixture():
    root = Path(__file__).resolve().parents[2]
    expected = json.loads((root / module.DATA).read_text())
    observed = {k: deepcopy(v) for k, v in expected.items() if k not in module.METADATA}
    return expected, observed


def test_exact_match():
    module.compare(*fixture())


@pytest.mark.parametrize("problem", ["interval", "version", "family", "seed", "extra", "publication"])
def test_scientific_changes_rejected(problem):
    expected, observed = fixture()
    if problem == "interval":
        observed["comparisons"][0]["metrics"]["f_score"]["paired_percentile_ci"][0] += .01
    elif problem == "version":
        observed["numpy_version"] = "other"
    elif problem == "family":
        observed["families"].pop()
    elif problem == "seed":
        observed["seed"] += 1
    elif problem == "extra":
        expected["unexplained"] = True
    else:
        expected["publication_ready"] = True
    with pytest.raises(ValueError):
        module.compare(expected, observed)


def test_existing_export_rejected(tmp_path):
    with pytest.raises(FileExistsError):
        module.reproduce(tmp_path, "HEAD", Path("python"), tmp_path)
