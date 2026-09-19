import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools import reproduce_qfo_sequence as module
from benchmark_tools.bootstrap_qfo_sequence import bootstrap

ROOT = Path(__file__).resolve().parents[2]


def evidence():
    directory = ROOT / "benchmark_tools/results"
    expected = json.loads((directory / module.RESULT).read_text())
    observed = bootstrap(json.loads((directory / module.COUNTS).read_text()))
    return expected, observed


def test_all_numerical_fields_match_real_admitted_evidence():
    module.compare(*evidence())


@pytest.mark.parametrize("field", sorted(module.NUMERICAL))
def test_no_numerical_field_is_ignored(field):
    expected, observed = evidence()
    observed[field] = None
    with pytest.raises(ValueError):
        module.compare(expected, observed)


@pytest.mark.parametrize("field,value", [("status", "admitted"), ("publication_ready", True),
    ("uncertainty_admitted", True), ("extra", 1)])
def test_numerical_reproduction_cannot_promote_admission(field, value):
    expected, observed = evidence()
    observed[field] = value
    with pytest.raises(ValueError):
        module.compare(expected, observed)


def test_worker_uses_local_pinned_counts_and_checks_inputs(tmp_path):
    export = tmp_path / "export"
    results = export / "benchmark_tools/results"
    results.mkdir(parents=True)
    for name in module.HASHES:
        shutil.copyfile(ROOT / "benchmark_tools/results" / name, results / name)
    module.worker(export, tmp_path / "result")
    generated = json.loads((tmp_path / "result/statistics.json").read_text())
    assert generated["replicates"] == 100000
    assert generated["uncertainty_admitted"] is False
    (results / module.COUNTS).write_text("{}")
    with pytest.raises(ValueError, match="Changed committed"):
        module.worker(export, tmp_path / "changed")
    assert not (tmp_path / "changed").exists()


def test_existing_export_refused(tmp_path):
    with pytest.raises(FileExistsError):
        module.reproduce(ROOT, "HEAD", Path("python"), tmp_path)
