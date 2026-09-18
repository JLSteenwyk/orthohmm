import copy
import json

import pytest

from benchmark_tools import export_qfo_sequence_scores as module


def report():
    participant = "ohmm_qfo_corrected_sequence_all_hits"
    conversion = dict(variant="all_hits", participant=participant,
        status="corrected_sequence_group_pairs_prepared_unscored",
        semantics="cross-species group-derived clique pairs",
        total_pairs=2, retained_pairs=2, expected_pairs=2, removed_mapping_pairs=0)
    return dict(status="corrected_sequence_assessment_admitted", variant="all_hits",
                accuracy_admitted=True, publication_ready=False, conversion=conversion,
                assessment=dict(participant=participant,
                    endpoints={key: dict(score=.5) for key in module.ENDPOINTS},
                    secondary_six_metric_mean=.5))


@pytest.mark.parametrize("problem", [None, "zero", "status", "variant", "participant", "semantics",
                                    "count", "boolean", "removed", "admission"])
def test_binding(problem):
    value = report()
    if problem == "zero":
        value["conversion"].update(total_pairs=0, retained_pairs=0, expected_pairs=0)
    elif problem == "status":
        value["status"] = "process_succeeded_pending_independent_admission"
    elif problem == "variant":
        value["variant"] = "unknown"
    elif problem == "participant":
        value["assessment"]["participant"] = "other"
    elif problem == "semantics":
        value["conversion"]["semantics"] = "native phylogenetic pairs"
    elif problem == "count":
        value["conversion"]["expected_pairs"] = 3
    elif problem == "boolean":
        value["conversion"]["removed_mapping_pairs"] = False
    elif problem == "removed":
        value["conversion"]["removed_mapping_pairs"] = 1
    elif problem == "admission":
        value["accuracy_admitted"] = False
    if problem in (None, "zero"):
        assert module.binding(value) == ("all_hits", "ohmm_qfo_corrected_sequence_all_hits")
    else:
        with pytest.raises(ValueError):
            module.binding(value)


def prepare(tmp_path, monkeypatch):
    value = report()
    def save(name, content):
        path = tmp_path / name
        path.write_text(json.dumps(content))
        return module.record(path)
    value["pairs_manifest"] = save("pairs.json", value["conversion"])
    value["execution_report"] = save("execution.json", dict(results=str(tmp_path / "results")))
    value["environment_manifest"] = save("env.json", dict(pipeline=str(tmp_path)))
    value["native_trace"] = save("trace.json", {})
    value["metric_files"] = []
    identity = save("admission.json", value)
    assessment = copy.deepcopy(value["assessment"])
    monkeypatch.setattr(module, "validate_directory", lambda *args: (assessment, []))
    return (identity["path"], identity["sha256"]), assessment


def test_export_preserves_pending_and_semantics(tmp_path, monkeypatch):
    source, _ = prepare(tmp_path, monkeypatch)
    out = tmp_path / "export"
    result = module.export([source], out)
    assert result["rows"][0]["status"] == "admitted"
    assert result["rows"][1]["status"] == "not_admitted"
    assert result["rows"][1]["scores"]["SwissTrees"] is None
    assert "pending" in (out / "scores.md").read_text()
    assert b"\r" not in (out / "scores.tsv").read_bytes()
    assert result["publication_ready"] is False
    with pytest.raises(FileExistsError):
        module.export([source], out)


def test_native_replay_mismatch_rejected(tmp_path, monkeypatch):
    source, assessment = prepare(tmp_path, monkeypatch)
    assessment["secondary_six_metric_mean"] = .6
    with pytest.raises(ValueError, match="Native metrics differ"):
        module.export([source], tmp_path / "export")
    assert not (tmp_path / "export").exists()


def test_duplicate_variant_rejected(tmp_path, monkeypatch):
    source, _ = prepare(tmp_path, monkeypatch)
    with pytest.raises(ValueError, match="Duplicate"):
        module.export([source, source], tmp_path / "export")


def test_retained_all_hit_scores_match_admitted_native_metrics():
    from pathlib import Path

    base = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    admission = module.read_frozen(base / "qfo_sequence_all_hits_assessment_admission_21826.json",
        "dc29e5410a274b5e633674bdf42c6578893bbc47912adc7b23ca6893ed917cd2")
    result = json.loads((base / "qfo_sequence_scores_20260918_v1/manifest.json").read_text())
    row, pending = result["rows"]
    assert module.binding(admission)[0] == row["variant"] == "all_hits"
    assert row["scores"] == {key: admission["assessment"]["endpoints"][key]["score"]
                             for key in module.ENDPOINTS}
    assert row["submitted_pairs"] == 11300151
    assert row["secondary_mean"] == pytest.approx(.6459776807897212)
    assert pending["status"] == "not_admitted"
    assert all(value is None for value in pending["scores"].values())
    module.check(result["source"])
    for item in result["outputs"]:
        module.check(item)


def test_both_retained_sequence_arms_match_admitted_endpoints():
    from pathlib import Path

    base = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    result = json.loads((base / "qfo_sequence_scores_20260918_v2/manifest.json").read_text())
    identities = [
        ("qfo_sequence_all_hits_assessment_admission_21826.json", "dc29e5410a274b5e633674bdf42c6578893bbc47912adc7b23ca6893ed917cd2", 11300151),
        ("qfo_sequence_top100_assessment_admission_21830.json", "439bfb1499e6c1f5e357ec27dc0fd02279d920a6a1a77f9483d008da4849c911", 11285357),
    ]
    for row, (name, sha, pairs) in zip(result["rows"], identities):
        admission = module.read_frozen(base / name, sha)
        assert module.binding(admission)[0] == row["variant"]
        assert row["status"] == "admitted"
        assert row["scores"] == {key: admission["assessment"]["endpoints"][key]["score"] for key in module.ENDPOINTS}
        assert row["secondary_mean"] == admission["assessment"]["secondary_six_metric_mean"]
        assert row["submitted_pairs"] == pairs
        assert row["removed_mapping_pairs"] == 0
        assert row["endpoint_details"] == admission["assessment"]["endpoints"]
    assert len(result["rows"]) == 2
    assert result["publication_ready"] is False
    module.check(result["source"])
    for item in result["outputs"]:
        module.check(item)
