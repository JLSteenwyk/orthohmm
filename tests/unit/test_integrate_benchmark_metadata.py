"""Preserve all selected results while integrating historical metadata scopes."""

from copy import deepcopy
import csv
import json
from pathlib import Path

import pytest

from benchmark_tools import integrate_benchmark_metadata as integration


def fixture_documents():
    base = {"status": "all_benchmark_provenance_register_partial", "complete_transitive_provenance": False,
            "scores_recomputed": False, "rows": []}
    supplement = {"status": "historical_three_kingdoms_metadata_supplement",
                  "historical_consumption_proven": False, "native_runs_repeated": False,
                  "scores_recomputed": False, "rows": []}
    scores = {"rows": [{"key": key, "scores": {"OrthoBench": .5, "ThreeKingdoms": .6,
                                             "QfO_secondary_mean": .7, "GO": .7}}
                       for key in integration.KEYS]}
    for dataset in integration.DATASETS:
        for key in integration.KEYS:
            pin = {"path": "/predictions/" + key, "bytes": 12, "sha256": "a" * 64}
            resources = [] if dataset == "QfO" and key.startswith("orthohmm") else [
                {"scope": "old recorded interval", "measurement": {"elapsed_seconds": 10},
                 "memory": {"value": None, "unit": "bytes", "scope": "unknown"}}]
            row = {"dataset": dataset, "key": key, "label": key, "declared_version": "retained",
                   "output_semantics": "groups", "scores": {"F1": .5 if dataset == "OrthoBench" else .6},
                   "commands": [], "resources": resources, "output_records": [pin], "gaps": ["old gap"]}
            if dataset == "QfO":
                row.update(scores={"GO": .7}, secondary_mean=.7)
            base["rows"].append(row)
            if dataset != "ThreeKingdoms" or key == "sonicparanoid_2_0_9":
                continue
            evidence = {"path": "/metrics.json", "bytes": 3, "sha256": "b" * 64}
            metrics = None
            times = []
            if key.startswith("orthohmm"):
                metrics = {"matched_output_manifest": [pin], "entrypoint_argv": ["python", "main.py"],
                           "harness_argv": ["python", "-m", "engine"], "measurement": {
                               "elapsed_seconds": 11.5, "user_seconds": 12, "system_seconds": 1},
                           "memory": {"value": 1234, "unit": "bytes", "scope": "sampled tree RSS"},
                           "cpu_scope": "recorded historical CPU"}
                files = {"metrics.json": evidence}
            else:
                names = (["bpo_converter.time.log", "downstream_time.log"] if key == "orthomcl_1_4" else
                         ["conversion_time.log"] if key.endswith("sequence_only") else ["time.log"])
                files = {}
                for name in names:
                    ref = {"path": "/historical/" + name, "bytes": 4, "sha256": "c" * 64}
                    files[name] = ref
                    times.append({"scope": "retained scope", "argv": ["tool", "--run"], "evidence": ref,
                                  "memory_scope": "GNU-time process RSS", "measurement": {
                                      "elapsed_seconds": 11.5, "user_seconds": 12, "system_seconds": 1,
                                      "max_process_rss_kib": 2345, "exit_status": 0}})
            supplement["rows"].append({"key": key, "declared_version": "retained", "selected_prediction": pin,
                "metrics": metrics, "timing_records": times, "evidence": files,
                "metadata_events": [{"key": "recovery", "value": "1"}, {"key": "recovery", "value": "2"}],
                "empty_timing_logs": ["time.log"] if key == "orthomcl_1_4" else [],
                "retained_text": {"source_commit.txt": "launch-not-completion"}})
    return base, supplement, scores


def test_every_original_field_and_input_document_preserved():
    documents = fixture_documents()
    before = deepcopy(documents)
    rows = integration.integrate(*documents)
    assert documents == before
    assert len(rows) == 24
    assert [{k: v for k, v in r.items() if k not in integration.NEW_FIELDS} for r in rows] == documents[0]["rows"]
    assert sum("historical_metadata_supplement" in r for r in rows) == 7
    sonic = next(r for r in rows if r["dataset"] == "ThreeKingdoms" and r["key"] == "sonicparanoid_2_0_9")
    assert not any(k in sonic for k in integration.NEW_FIELDS)


def test_resource_placeholders_scopes_and_recovery_events_preserved():
    rows = integration.integrate(*fixture_documents())
    resources = integration.resource_rows(rows)
    assert len(resources) == 32  # 24 base rows, including two placeholders, plus eight supplemental intervals.
    missing = [r for r in resources if r["wall_measurement_status"] == "unavailable"]
    assert len(missing) == 2 and all(r["dataset"] == "QfO" and r["wall_seconds"] is None for r in missing)
    assert all(r["independent_repeat"] is False for r in resources)
    mcl = next(r for r in rows if r["dataset"] == "ThreeKingdoms" and r["key"] == "orthomcl_1_4")
    assert [r["value"] for r in mcl["historical_metadata_supplement"]["metadata_events"]] == ["1", "2"]
    assert "excludes BLAST/BPO" in mcl["supplemental_resources"][1]["scope"]
    checkpoint = next(r for r in rows if r["dataset"] == "ThreeKingdoms" and r["key"].endswith("sequence_only"))
    assert "not separate sequence-only inference" in checkpoint["supplemental_resources"][0]["scope"]
    fastoma = next(r for r in rows if r["dataset"] == "ThreeKingdoms" and r["key"] == "fastoma_0_3_5")
    assert "exclude aggregate Docker tasks" in fastoma["supplemental_resources"][0]["scope"]


@pytest.mark.parametrize("damage", ["duplicate", "missing", "score", "secondary", "prediction", "version",
                                    "sonic", "scope", "argv", "failure", "timing_pin", "empty_output", "already_added"])
def test_wrong_bindings_and_scopes_rejected(damage):
    base, supplement, scores = fixture_documents()
    if damage == "duplicate":
        base["rows"][-1] = deepcopy(base["rows"][-2])
    elif damage == "missing":
        base["rows"].pop()
    elif damage == "score":
        base["rows"][0]["scores"]["F1"] = .51
    elif damage == "secondary":
        base["rows"][8]["secondary_mean"] = .71
    elif damage == "prediction":
        supplement["rows"][0]["selected_prediction"] = {"path": "/wrong", "bytes": 2, "sha256": "d" * 64}
    elif damage == "version":
        supplement["rows"][0]["declared_version"] = "different"
    elif damage == "sonic":
        supplement["rows"][0]["key"] = "sonicparanoid_2_0_9"
    elif damage == "scope":
        supplement["historical_consumption_proven"] = True
    elif damage == "argv":
        supplement["rows"][0]["metrics"]["entrypoint_argv"] = "python main.py"
    elif damage == "failure":
        supplement["rows"][2]["timing_records"][0]["measurement"]["exit_status"] = 1
    elif damage == "timing_pin":
        supplement["rows"][2]["evidence"]["time.log"] = {"path": "/wrong"}
    elif damage == "empty_output":
        base["rows"][16]["output_records"] = []
    else:
        base["rows"][0]["supplemental_commands"] = []
    with pytest.raises(ValueError):
        integration.integrate(base, supplement, scores)


@pytest.mark.parametrize("value", [None, True, -1, float("nan"), float("inf")])
def test_bad_resource_values_rejected(value):
    documents = fixture_documents()
    documents[1]["rows"][0]["metrics"]["measurement"]["elapsed_seconds"] = value
    with pytest.raises(ValueError):
        integration.integrate(*documents)


def test_existing_output_untouched(tmp_path):
    output = tmp_path / "result"
    output.mkdir()
    with pytest.raises(FileExistsError):
        integration.write(tmp_path, output)
    assert list(output.iterdir()) == []


def test_actual_24_row_readback_and_every_tsv_field():
    repo = Path(__file__).resolve().parents[2]
    output = repo / "benchmark_tools/results/all_benchmark_metadata_integrated_20261004_v3"
    if not output.exists():
        pytest.skip("Actual integrated report not yet available")
    report = json.loads((output / "register.json").read_text())
    current_source = integration.record(repo / "benchmark_tools/integrate_benchmark_metadata.py")
    assert all(current_source[k] == report["source"][k] for k in ("bytes", "sha256"))
    documents = {key: json.loads((repo / "benchmark_tools/results" / spec[0]).read_text())
                 for key, spec in integration.SOURCES.items()}
    for key, spec in integration.SOURCES.items():
        observed = integration.record(repo / "benchmark_tools/results" / spec[0])
        assert observed["sha256"] == spec[1]
        inherited = next(r for r in report["inputs"] if r["sha256"] == spec[1])
        assert observed["bytes"] == inherited["bytes"]
    assert report["rows"] == integration.integrate(documents["register"], documents["supplement"], documents["scores"])
    assert report["resource_intervals"] == integration.resource_rows(report["rows"])
    assert len(report["resource_intervals"]) == 37
    assert sum(r["wall_seconds"] is not None for r in report["resource_intervals"]) == 34
    assert (output / "register.md").read_text() == integration.render(report)
    with (output / "resources.tsv").open(newline="") as stream:
        records = list(csv.DictReader(stream, delimiter="\t"))
    expected = [{key: "NA" if value is None else str(value) for key, value in row.items()}
                for row in report["resource_intervals"]]
    assert records == expected
    assert report["previous_register_fields_unchanged"] is True
    for field in ("historical_input_consumption_proven", "complete_transitive_provenance",
                  "controlled_comparative_resources", "publication_ready", "native_inference_or_scoring_repeated"):
        assert report[field] is False
