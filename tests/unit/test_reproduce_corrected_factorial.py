from copy import deepcopy
import json
import math
from pathlib import Path
import shutil

import pytest

from benchmark_tools import reproduce_corrected_factorial as module
from benchmark_tools.bootstrap_qfo_factorial import bootstrap
from benchmark_tools import bootstrap_qfo_factorial as engine

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def evidence():
    directory = ROOT / "benchmark_tools/results"
    expected = json.loads((directory / module.RESULT).read_text())
    counts = json.loads((directory / module.COUNTS).read_text())
    observed = bootstrap({**counts, "status": "qfo_factorial_swiss_counts_verified"})
    return expected, observed


@pytest.mark.retained_benchmark
def test_all_numerical_fields_match(evidence):
    module.compare(*evidence)


def canonical(expected):
    observed = {key: deepcopy(expected[key]) for key in module.NUMERICAL | module.COMMON}
    observed["status"] = "paired_qfo_factorial_swiss_intervals"
    return observed


def roundoff_differences(left, right, path=""):
    """Diagnose only test-platform roundoff; never admit a scientific replay."""
    assert type(left) is type(right), path
    differences = []
    if isinstance(left, dict):
        assert left.keys() == right.keys(), path
        for key in left:
            differences.extend(roundoff_differences(left[key], right[key], path + "/" + str(key)))
    elif isinstance(left, list):
        assert len(left) == len(right), path
        for index, (a, b) in enumerate(zip(left, right)):
            differences.extend(roundoff_differences(a, b, path + "/" + str(index)))
    elif isinstance(left, float):
        assert math.isfinite(left) and math.isfinite(right), path
        error = abs(left - right)
        assert error <= 64 * math.ulp(1.0), (path, left, right)
        if error:
            differences.append((path, error))
    else:
        assert left == right, path
    return differences


def test_platform_roundoff_is_bounded_but_never_promoted_to_exact_replay(evidence, record_testsuite_property):
    expected, observed = evidence
    left, right = ({key: value[key] for key in module.NUMERICAL} for value in evidence)
    differences = roundoff_differences(left, right)
    record_testsuite_property("corrected_factorial_floating_differences", len(differences))
    record_testsuite_property("corrected_factorial_maximum_absolute_roundoff", max(
        (error for _, error in differences), default=0))
    if differences:
        with pytest.raises(ValueError, match="scientific statistics differ"):
            module.compare(expected, observed)
    else:
        module.compare(expected, observed)


@pytest.mark.parametrize("left,right", [(1, 2), (True, False), (1., 1.000001),
    (1., float("nan")), ({"a": 1}, {"b": 1}), ([1], [1, 2]), ("2.2.6", "2.2.7")])
def test_roundoff_diagnostic_rejects_changed_metadata_counts_and_statistics(left, right):
    with pytest.raises(AssertionError):
        roundoff_differences(left, right)


def test_strict_comparison_rejects_even_one_ulp(evidence):
    expected, _ = evidence
    observed = canonical(expected)
    module.compare(expected, observed)
    point = next(iter(observed["point_estimates"].values()))
    point["F1"] = math.nextafter(point["F1"], math.inf)
    with pytest.raises(ValueError, match="scientific statistics differ"):
        module.compare(expected, observed)


@pytest.mark.parametrize("field", sorted(module.NUMERICAL))
def test_no_numerical_field_ignored(evidence, field):
    expected, _ = evidence
    observed = canonical(expected)
    module.compare(expected, observed)
    observed[field] = None
    with pytest.raises(ValueError):
        module.compare(expected, observed)


@pytest.mark.parametrize("field,value", [("status", "admitted"), ("publication_ready", True), ("extra", 1)])
def test_no_admission_promotion(evidence, field, value):
    expected, _ = evidence
    observed = canonical(expected)
    module.compare(expected, observed)
    observed[field] = value
    with pytest.raises(ValueError):
        module.compare(expected, observed)


def test_wrong_release_rejected(evidence):
    expected, _ = deepcopy(evidence)
    observed = canonical(expected)
    expected["input_release"] = "original"
    with pytest.raises(ValueError):
        module.compare(expected, observed)


@pytest.fixture
def export(tmp_path):
    export = tmp_path / "export"
    results = export / "benchmark_tools/results"
    for name in module.HASHES:
        path = results / name
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / "benchmark_tools/results" / name, path)
    engine = export / "benchmark_tools/bootstrap_qfo_factorial.py"
    shutil.copyfile(ROOT / "benchmark_tools/bootstrap_qfo_factorial.py", engine)
    return export


def test_worker_pins_inputs_and_engine_with_canonical_fixture(tmp_path, export, evidence, monkeypatch):
    expected, _ = evidence
    calls = []

    def canonical_bootstrap(report, **kwargs):
        calls.append((deepcopy(report), kwargs))
        return canonical(expected)

    monkeypatch.setattr(engine, "bootstrap", canonical_bootstrap)
    module.worker(export, tmp_path / "output")
    assert json.loads((tmp_path / "output/statistics.json").read_text())["multiplicity_endpoints"] == 42
    assert len(calls) == 1 and calls[0][1] == dict(replicates=100000, seed=20260922)
    counts = json.loads((export / "benchmark_tools/results" / module.COUNTS).read_text())
    assert calls[0][0] == dict(counts, status="qfo_factorial_swiss_counts_verified")
    with pytest.raises(FileExistsError):
        module.worker(export, tmp_path / "output")
    (export / "benchmark_tools/bootstrap_qfo_factorial.py").write_text("changed")
    with pytest.raises(ValueError, match="engine"):
        module.worker(export, tmp_path / "changed")
    assert not (tmp_path / "changed").exists()


def test_worker_refuses_roundoff_without_writing_output(tmp_path, export, evidence, monkeypatch):
    expected, _ = evidence
    observed = canonical(expected)
    point = next(iter(observed["point_estimates"].values()))
    point["F1"] = math.nextafter(point["F1"], math.inf)
    monkeypatch.setattr(engine, "bootstrap", lambda *a, **k: deepcopy(observed))
    with pytest.raises(ValueError, match="scientific statistics differ"):
        module.worker(export, tmp_path / "output")
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("name", sorted(module.HASHES))
def test_worker_refuses_changed_input_bytes_before_bootstrap(tmp_path, export, name, monkeypatch):
    (export / "benchmark_tools/results" / name).write_text("changed input\n")
    monkeypatch.setattr(engine, "bootstrap", lambda *a, **k: pytest.fail("must validate inputs first"))
    with pytest.raises(ValueError, match="Changed committed analysis input"):
        module.worker(export, tmp_path / "output")
    assert not (tmp_path / "output").exists()


def test_existing_export_refused(tmp_path):
    with pytest.raises(FileExistsError):
        module.reproduce(ROOT, "HEAD", Path("python"), tmp_path)
