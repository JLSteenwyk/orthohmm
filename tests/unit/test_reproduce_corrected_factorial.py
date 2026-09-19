from copy import deepcopy
import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools import reproduce_corrected_factorial as module
from benchmark_tools.bootstrap_qfo_factorial import bootstrap

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def evidence():
    directory = ROOT / "benchmark_tools/results"
    expected = json.loads((directory / module.RESULT).read_text())
    counts = json.loads((directory / module.COUNTS).read_text())
    observed = bootstrap({**counts, "status": "qfo_factorial_swiss_counts_verified"})
    return expected, observed


def test_all_numerical_fields_match(evidence):
    module.compare(*evidence)


@pytest.mark.parametrize("field", sorted(module.NUMERICAL))
def test_no_numerical_field_ignored(evidence, field):
    expected, observed = deepcopy(evidence)
    observed[field] = None
    with pytest.raises(ValueError):
        module.compare(expected, observed)


@pytest.mark.parametrize("field,value", [("status", "admitted"), ("publication_ready", True), ("extra", 1)])
def test_no_admission_promotion(evidence, field, value):
    expected, observed = deepcopy(evidence)
    observed[field] = value
    with pytest.raises(ValueError):
        module.compare(expected, observed)


def test_wrong_release_rejected(evidence):
    expected, observed = deepcopy(evidence)
    expected["input_release"] = "original"
    with pytest.raises(ValueError):
        module.compare(expected, observed)


def test_worker_pins_inputs_and_engine(tmp_path):
    export = tmp_path / "export"
    results = export / "benchmark_tools/results"
    for name in module.HASHES:
        path = results / name
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / "benchmark_tools/results" / name, path)
    engine = export / "benchmark_tools/bootstrap_qfo_factorial.py"
    shutil.copyfile(ROOT / "benchmark_tools/bootstrap_qfo_factorial.py", engine)
    module.worker(export, tmp_path / "output")
    assert json.loads((tmp_path / "output/statistics.json").read_text())["multiplicity_endpoints"] == 42
    with pytest.raises(FileExistsError):
        module.worker(export, tmp_path / "output")
    engine.write_text("changed")
    with pytest.raises(ValueError, match="engine"):
        module.worker(export, tmp_path / "changed")
    assert not (tmp_path / "changed").exists()


def test_existing_export_refused(tmp_path):
    with pytest.raises(FileExistsError):
        module.reproduce(ROOT, "HEAD", Path("python"), tmp_path)
