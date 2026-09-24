import json

import pytest

from benchmark_tools import export_qfo_complete_comparison as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_export_qfo_corrected_factorial import fixture
from tests.unit.test_export_qfo_corrected_comparison import fixture as comparator_fixture, orthomcl_fixture, recovered_fixture


@pytest.mark.parametrize("method", ["proteinortho", "sonic", "orthomcl"])
def test_existing_comparators_keep_exact_admitted_fields(method):
    report, conversion = orthomcl_fixture() if method == "orthomcl" else comparator_fixture(method)
    assert module.extract(report, conversion) == module.comparator(report, conversion)


@pytest.mark.parametrize("index", range(8))
def test_only_frozen_publication_cells(index):
    report = fixture(index)
    if index not in (4, 7):
        with pytest.raises(ValueError, match="not a frozen publication"):
            module.extract(report, report["conversion"])
        return
    row = module.extract(report, report["conversion"])
    assert row["key"] == module.PUBLICATION_CELLS[index]
    assert row["scores"]["SwissTrees"] == .5
    assert row["submitted_pairs"] == row["retained_pairs"] == 2
    assert row["prediction_semantics"] == report["conversion"]["semantics"]


def inputs(tmp_path, monkeypatch):
    replay = tmp_path / "replay.json"
    replay.write_text(json.dumps({"status": "corrected_checked_replay_admitted"}))
    monkeypatch.setattr(module, "REPLAY_SHA", record(replay)["sha256"])
    report = fixture(4)
    pairs = tmp_path / "pairs.json"
    pairs.write_text(json.dumps(report["conversion"]))
    report["pairs_manifest"] = record(pairs)
    admission = tmp_path / "admission.json"
    admission.write_text(json.dumps(report))
    return replay, admission, pairs, report


def test_partial_export(tmp_path, monkeypatch):
    replay, admission, pairs, report = inputs(tmp_path, monkeypatch)
    output = tmp_path / "table"
    sources = [(admission, record(admission)["sha256"])]
    result = module.export(sources, replay, output)
    assert result["admitted_methods"] == 1 and len(result["methods"]) == 8
    assert result["methods"][0]["factorial_cell"] == "p1_c0_r0"
    assert all(r["secondary_mean"] is None for r in result["methods"][1:])
    assert "cross-species group-derived clique pairs" in (output / "scores.md").read_text()
    assert "Mapping losses" in (output / "scores.tsv").read_text()
    with pytest.raises(FileExistsError):
        module.export(sources, replay, output)


def test_recovered_export_joins_without_relabeling(tmp_path, monkeypatch):
    replay, existing, _, _ = inputs(tmp_path, monkeypatch)
    report, conversion = recovered_fixture()
    pairs = tmp_path / "recovered_pairs.json"
    pairs.write_text(json.dumps(conversion))
    report["pairs_manifest"] = record(pairs)
    path = tmp_path / "recovered_admission.json"
    path.write_text(json.dumps(report))
    output = tmp_path / "table"
    result = module.export([(p, record(p)["sha256"]) for p in (existing, path)], replay, output)
    assert result["admitted_methods"] == 2
    recovered = next(r for r in result["methods"] if r["key"] == "orthomcl_1_4")
    assert recovered["participant"] == "qfo_corrected_orthomcl_recovered"
    assert recovered["query_coverage"]["failed_queries"] == ["retained"]
    assert module.RECOVERY_NOTE in (output / "scores.md").read_text()


@pytest.mark.parametrize("problem", ["replay", "report", "conversion", "duplicate", "unadmitted", "wrong_cell"])
def test_fail_closed(tmp_path, monkeypatch, problem):
    replay, admission, pairs, report = inputs(tmp_path, monkeypatch)
    sources = [(admission, record(admission)["sha256"])]
    if problem == "replay":
        replay.write_text("{}")
    elif problem == "report":
        admission.write_text("{}")
    elif problem == "conversion":
        pairs.write_text("{}")
    elif problem == "duplicate":
        sources *= 2
    else:
        if problem == "unadmitted":
            report["accuracy_admitted"] = False
        else:
            report["cell"] = "p0_c0_r0"
        admission.write_text(json.dumps(report))
        sources = [(admission, record(admission)["sha256"])]
    output = tmp_path / "out"
    with pytest.raises(ValueError):
        module.export(sources, replay, output)
    assert not output.exists()


def test_broken_output_symlink(tmp_path):
    output = tmp_path / "out"
    output.symlink_to(tmp_path / "absent")
    with pytest.raises(FileExistsError):
        module.export([], tmp_path / "missing", output)
