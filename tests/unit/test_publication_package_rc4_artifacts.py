"""Read back actual delivered rc4 anchors and the copied reporting execution."""

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"


def read(name):
    return json.loads((RESULTS / name).read_bytes())


def test_selected_candidate_and_actual_index_match_exact_anchors():
    expected = {
        "publication_package_rc4_selection_20261004.json": "77784822a38c59257050752a02aa87dc3df48a6c9a3f07a764e06dd4df097283",
        "publication_package_rc4_index_20261004.json": "550d1b06875a26d5719376d5f8562055f68b87ba5d05619b1273146c4c51cf71",
    }
    for path, sha in expected.items():
        assert hashlib.sha256((RESULTS / path).read_bytes()).hexdigest() == sha
    selection = read("publication_package_rc4_selection_20261004.json")
    index = read("publication_package_rc4_index_20261004.json")
    assert len(selection["files"]) == 181 and len(index["files"]) == 183
    assert index["workflow_revision"] == selection["workflow_revision"] == "cd1894ca3979ec10f67ba5d71e5159908c8a61b3"
    assert index["scientific_revision"] == selection["scientific_revision"] == "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
    assert index["publication_ready"] is index["public_archive_uploaded"] is False


def test_all_139_prior_indexed_payload_identities_survive_with_explicit_relocations():
    previous = read("publication_package_rc3_index_20261004.json")
    current = {r["path"]: r for r in read("publication_package_rc4_index_20261004.json")["files"]}
    relocated = {"README.md": "evidence/PUBLICATION_PACKAGE_RC3_20261004.md",
                 "PACKAGE_SELECTION.json": "history/rc3/PACKAGE_SELECTION.json",
                 "bundle_publication_package.py": "history/rc3/bundle_publication_package.py"}
    assert len(previous["files"]) == 139
    for row in previous["files"]:
        target = relocated.get(row["path"], row["path"])
        assert {k: current[target][k] for k in ("bytes", "sha256")} == {k: row[k] for k in ("bytes", "sha256")}
    review = current["current-review/document.pdf"]
    assert review["sha256"] == "feb0b085def6eed4d2c111668e8ebb6571c8d55e38c7c3a80cdc01b34f567a9b"
    assert current["review/document.pdf"]["sha256"] == next(
        r["sha256"] for r in previous["files"] if r["path"] == "review/document.pdf")


def test_copied_addenda_reporting_outputs_preserve_quantitative_scope():
    report = read("publication_addenda_reporting_replay_20261004.json")
    assert hashlib.sha256((RESULTS / "publication_addenda_reporting_replay_20261004.json").read_bytes()).hexdigest() == "bc0a24dc4aba948a7f0b830c66dc428ff3e7537224e7e58ce1a9160116b9fa63"
    assert report["raw_or_historical_paths_accessed"] is report["native_or_scoring_repeated"] is report["publication_ready"] is False
    summary = report["summary"]
    assert summary["exposure"]["files"] == 1899
    assert summary["exposure"]["families"] == 88 and summary["exposure"]["associations"] == 9156
    assert summary["register"]["metric_positions"] == 64 and summary["register"]["secondary_mean_positions"] == 8
    assert summary["register"]["stage_associations"] == 5 and summary["register"]["distinct_stage_observations"] == 4
    assert summary["register"]["original_full_cost_cells_unavailable"] == 16
    assert all(r["unattached_count"] == 1 for r in summary["candidate_trace"]["fixture_cases"])
