import json
from pathlib import Path

import pytest

from benchmark_tools.audit_publication_figures import inspect_manifest, records, PANELS
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture(tmp_path):
    directory = tmp_path / "figures"
    directory.mkdir()
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = directory / ("figure." + extension)
        path.write_bytes(b"synthetic format placeholder")
        outputs.append(record(path))
    source = tmp_path / "source.py"
    source.write_text("pass\n")
    data = {"outputs": outputs, "source": record(source)}
    path = directory / "manifest.json"
    path.write_text(json.dumps(data))
    return path, data


def test_checks_identity_and_separates_tracking(tmp_path):
    path, data = fixture(tmp_path)
    result = inspect_manifest(path, tmp_path, {"source.py"})
    assert result["status"] == "all_recorded_bytes_match"
    assert result["output_count"] == 3
    assert sum(row["tracked_in_main_repository"] for row in result["files"]) == 1
    assert len(result["files"]) == 4


@pytest.mark.parametrize("mutation", ["changed", "missing"])
def test_preserves_integrity_failures(tmp_path, mutation):
    path, data = fixture(tmp_path)
    source = Path(data["source"]["path"])
    if mutation == "changed":
        source.write_text("changed\n")
    else:
        source.unlink()
    result = inspect_manifest(path, tmp_path, set())
    assert result["status"] == "integrity_failure"
    assert result["files"][-1]["status"] == mutation


@pytest.mark.parametrize("mutation", ["duplicate", "format", "outside", "conflict", "schema"])
def test_rejects_invalid_inventory(tmp_path, mutation):
    path, data = fixture(tmp_path)
    if mutation == "duplicate":
        data["outputs"].append(data["outputs"][0])
    elif mutation == "format":
        data["outputs"].pop()
    elif mutation == "outside":
        data["outputs"][0]["path"] = str(tmp_path / "elsewhere.png")
    elif mutation == "conflict":
        data["extra"] = {**data["source"], "bytes": 999}
    else:
        data["source"]["unexpected"] = True
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        inspect_manifest(path, tmp_path, set())


def test_nested_record_inventory():
    item = {"path": "/tmp/a", "bytes": 1, "sha256": "a"}
    assert list(records({"nested": [item, {"deeper": item}], "text": "ignored"})) == [item, item]


def test_original_qfo_factorial_is_explicitly_retained():
    assert "qfo_factorial_swiss_figure_20260918" in PANELS
    assert "figures_dgx_descriptive_20260918" in PANELS
    assert "figures_qfo_hit_coverage_20260918" in PANELS
    assert "figures_qfo_sequence_search_20260918" in PANELS
    assert len(set(PANELS)) == len(PANELS)
    assert not any("corrected_factorial" in panel for panel in PANELS)


@pytest.mark.parametrize("release", ["corrected", "original", None])
def test_corrected_scope_is_separate_and_release_checked(tmp_path, monkeypatch, release):
    from benchmark_tools import audit_publication_figures as module
    original, data = fixture(tmp_path)
    directory = tmp_path / "benchmark_tools/results" / module.CORRECTED_PANELS[0]
    directory.parent.mkdir(parents=True)
    original.parent.rename(directory)
    data["input_release"] = release
    data["outputs"] = [record(directory / Path(row["path"]).name) for row in data["outputs"]]
    (directory / "manifest.json").write_text(json.dumps(data))
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "")
    if release != "corrected":
        with pytest.raises(ValueError, match="corrected-release"):
            module.audit(tmp_path, "corrected-factorial")
    else:
        result = module.audit(tmp_path, "corrected-factorial")
        assert result["status"] == "retained_figure_bytes_verified"
        assert result["scope"] == "corrected-factorial"
        assert len(result["panels"]) == 1
        assert result["total_output_records"] == 3
        assert not result["publication_ready"]


def test_unknown_audit_scope_rejected(tmp_path):
    from benchmark_tools.audit_publication_figures import audit
    with pytest.raises(ValueError, match="scope"):
        audit(tmp_path, "combined-implicit")


def test_current_corrected_scope_does_not_change_old_inventories():
    from benchmark_tools.audit_publication_figures import CORRECTED_PANELS, CORRECTED_CURRENT_PANELS
    assert CORRECTED_PANELS == ("qfo_corrected_factorial_figures_20260919",)
    assert CORRECTED_CURRENT_PANELS == (*CORRECTED_PANELS, "corrected_swiss_strata_figure_21981",
                                       "corrected_swiss_comparison_figure_21987")
    assert not set(CORRECTED_CURRENT_PANELS).intersection(PANELS)


@pytest.mark.parametrize("kind", ["strata", "comparison", "fastoma", "complete"])
@pytest.mark.parametrize("problem", [None, "status", "admission", "ready", "endpoints", "duplicate"])
def test_new_corrected_panel_source_binding(tmp_path, kind, problem):
    from benchmark_tools.audit_publication_figures import validate_corrected_panel
    source = tmp_path / "qfo_corrected_comparator_uncertainty_21987.json"
    if kind == "fastoma":
        source = tmp_path / "qfo_fastoma_swiss_uncertainty_22098.json"
    elif kind == "complete":
        source = tmp_path / "qfo_recovered_swiss_uncertainty_22178.json"
    result = {"status": "corrected_swiss_primary_stratified_intervals" if kind == "strata"
              else "corrected_swiss_comparison_intervals_audited", "scientific_inputs_admitted": True,
              "uncertainty_admitted": True, "publication_ready": False}
    if kind == "complete":
        result.update(complete_panel=True, estimated_contrasts=8, multiplicity_endpoints=24,
                      point_estimates={str(i): {} for i in range(8)})
    if problem == "status":
        result["status"] = "paired_swiss_comparator_intervals"
    elif problem == "admission":
        result["scientific_inputs_admitted"] = False
    elif problem == "ready":
        result["publication_ready"] = True
    source.write_text(json.dumps(result))
    data = {"publication_ready": False, "endpoints": 27 if kind == "strata" else 24,
            "results": record(source), "inputs": [record(source)], "status": "corrected_swiss_comparison_rendered"}
    if problem == "endpoints":
        data["endpoints"] = 0
    elif problem == "duplicate":
        if kind == "strata":
            data["publication_ready"] = True
        else:
            data["inputs"] *= 2
    panel = "corrected_swiss_strata_figure_21981" if kind == "strata" else "corrected_swiss_comparison_figure_21987"
    if kind == "fastoma":
        panel = "corrected_swiss_comparison_figure_20260923"
    elif kind == "complete":
        panel = "corrected_swiss_comparison_figure_20260926"
    if problem:
        with pytest.raises(ValueError):
            validate_corrected_panel(panel, data)
    else:
        validate_corrected_panel(panel, data)


def test_fastoma_scope_rejects_old_comparator_source(tmp_path):
    from benchmark_tools.audit_publication_figures import validate_corrected_panel, CORRECTED_FASTOMA_PANELS
    assert CORRECTED_FASTOMA_PANELS == ("corrected_swiss_comparison_figure_20260923",)
    data = {"publication_ready": False, "endpoints": 24,
            "status": "corrected_swiss_comparison_rendered",
            "inputs": [{"path": str(tmp_path / "qfo_corrected_comparator_uncertainty_21987.json")}]}
    with pytest.raises(ValueError, match="one corrected comparator"):
        validate_corrected_panel(CORRECTED_FASTOMA_PANELS[0], data)


@pytest.mark.parametrize("field, value", [
    ("complete_panel", False), ("estimated_contrasts", 7),
    ("multiplicity_endpoints", 21), ("point_estimates", {}),
])
def test_complete_scope_rejects_partial_panel(tmp_path, field, value):
    from benchmark_tools.audit_publication_figures import validate_corrected_panel, CORRECTED_COMPLETE_PANELS
    source = tmp_path / "qfo_recovered_swiss_uncertainty_22178.json"
    result = dict(status="corrected_swiss_comparison_intervals_audited",
                  scientific_inputs_admitted=True, uncertainty_admitted=True, publication_ready=False,
                  complete_panel=True, estimated_contrasts=8, multiplicity_endpoints=24,
                  point_estimates={str(i): {} for i in range(8)})
    result[field] = value
    source.write_text(json.dumps(result))
    data = dict(publication_ready=False, endpoints=24, status="corrected_swiss_comparison_rendered",
                inputs=[record(source)])
    with pytest.raises(ValueError, match="complete eight-method"):
        validate_corrected_panel(CORRECTED_COMPLETE_PANELS[0], data)
