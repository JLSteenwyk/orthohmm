"""Retained selected artifacts and disclosure checks, not another scientific run."""

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"
PINS = {
    "swiss_ordered_pfam_features_20261007_v1.json": "c002bf9d61b523545c9f46002e147ae5c04bacc9335d48c04886ffe549f0810f",
    "swiss_ordered_pfam_feature_readback_20261007_v1.json": "275b89720a65d0c5d0a9afe70a806c53a3d9fbb8014d99d3d4ba361f78bf56bc",
    "native_swiss_ordered_pfam_20261007_v1/report.json": "f3389ad993f7ffbf5a5ee3fe2e2e8bac061120dc10947237a082d2c9e13ea9b0",
    "native_swiss_ordered_pfam_readback_20261007_v1.json": "db942129441adcc470164ddef99a8ed0c157f299edd0cd76b0879c2066a2caba",
    "native_swiss_ordered_pfam_20261007_v1/TABLE.md": "f220d342a4c90109b9569f4c7f311de00c7a60c083c7aa1434f5db9fb49c1829",
    "native_swiss_ordered_pfam_20261007_v1/scores.tsv": "ecf046dac73309cb2e34e945dc3cfeb0b53c727c3953d88d8b7002060309445d",
    "native_swiss_ordered_pfam_20261007_v1/differences.tsv": "45f4665443b8f130e4006c5df8ecca4f1f9a3300183a2c7cef8ee93d668d7c08",
}


def document(name):
    return json.loads((RESULTS / name).read_text())


def test_selected_bytes_bindings_complete_universe_and_scope():
    for name, sha in PINS.items():
        assert hashlib.sha256((RESULTS / name).read_bytes()).hexdigest() == sha
    features = document("swiss_ordered_pfam_features_20261007_v1.json")
    reader = document("swiss_ordered_pfam_feature_readback_20261007_v1.json")
    projected = document("native_swiss_ordered_pfam_20261007_v1/report.json")
    rational = document("native_swiss_ordered_pfam_readback_20261007_v1.json")
    assert reader["report"]["sha256"] == PINS["swiss_ordered_pfam_features_20261007_v1.json"]
    assert rational["report"]["sha256"] == PINS["native_swiss_ordered_pfam_20261007_v1/report.json"]
    assert len(features["genes"]) == 563 and len(features["families"]) == 18
    assert projected["memberships"] == features["memberships"]
    original_counts = json.loads(Path(projected["inputs"]["counts"]["path"]).read_text())
    assert projected["family_rows"] == original_counts["family_rows"] and len(projected["family_rows"]) == 54
    assert len(projected["rows"]) == 12 and len(projected["differences"]) == 8
    assert reader["status"] == "ordered_features_verified" and rational["status"] == "ordered_projection_verified"
    assert rational["human_numeric_rows_checked"] == 20
    assert projected["cells"][1]["timing_eligible"] is False and projected["cells"][1]["timing_admitted"] is False
    for report in (features, reader, projected, rational):
        for key in ("publication_ready", "independent_confirmation", "new_uncertainty", "scientific_timings_admitted"):
            assert report[key] is False
        assert report["new_bootstrap_draws"] == 0


def test_annotation_comparability_does_not_become_biological_rearrangement_truth():
    features = document("swiss_ordered_pfam_features_20261007_v1.json")
    reader = document("swiss_ordered_pfam_feature_readback_20261007_v1.json")
    assert reader["protein_states"] == dict(usable=434, order_ambiguous=124, zero_pfam=5)
    assert reader["same_multiset_comparable_pairs"] == 5564
    assert reader["same_multiset_order_discordant_pairs"] == 0
    assert len(features["bins"]["all_members_usable_same_signature"]) == 4
    assert features["bins"]["all_members_usable_multiple_signatures"] == ["Clusterin"]
    assert len(features["bins"]["some_members_unusable"]) == 13
    signatures = {tuple(features["genes"][g]["ordered_signature"]) for g in features["families"]["Clusterin"]["members"]}
    assert signatures == {("pfam_Clusterin",), ("pfam_Clusterin", "pfam_Clusterin")}
    text = " ".join((RESULTS / "SWISS_ORDERED_PFAM_RESULT_20261007.md").read_text().split())
    for phrase in ("not a rearrangement", "not independent observations", "not biological loss", "single family",
                   "shared-host observations", "No negative finding proves", "not literal full"):
        assert phrase in text


def test_actual_commands_are_four_once_only_isolated_stages_with_diagnostic_times():
    execution = document("swiss_ordered_pfam_execution_20261007_v1.json")
    assert execution["source_commit"] == "d09e2bc41e12020c6217e28e22811c67477eaee4"
    assert execution["counts_entered_only_after_full_feature_readback"] is True
    assert execution["runtime"]["python"] == "3.10.13" and execution["runtime"]["biopython"] == "1.87"
    assert [r["stage"] for r in execution["selected_invocations"]] == ["features", "feature_readback", "projection", "projection_readback"]
    for row in execution["selected_invocations"]:
        assert row["invocations"] == 1 and row["exit_code"] == 0
        assert " -I -B " in row["command"] and "OPENBLAS_NUM_THREADS=1" in row["command"]
        assert execution["source_commit"] in row["command"]
        assert (ROOT / row["time_file"]).read_text().startswith("exit=0\nelapsed_seconds=")
    for key, name in (("prepare_swiss_ordered_pfam", "prepare_swiss_ordered_pfam.py"),
                      ("readback_swiss_ordered_pfam", "readback_swiss_ordered_pfam.py")):
        assert hashlib.sha256((ROOT / "benchmark_tools" / name).read_bytes()).hexdigest() == execution["source_hashes"][key]
    assert execution["scientific_timings_admitted"] is False and execution["publication_ready"] is False


def test_result_summary_matches_generated_conditional_f1_points():
    projected = document("native_swiss_ordered_pfam_20261007_v1/report.json")
    text = (RESULTS / "SWISS_ORDERED_PFAM_RESULT_20261007.md").read_text()
    labels = dict(all="All", all_members_usable_same_signature="All usable, same signature",
                  all_members_usable_multiple_signatures="All usable, multiple signatures", some_members_unusable="Some members unusable")
    differences = {(r["stratum"], r["contrast"]): r for r in projected["differences"]}
    for stratum, label in labels.items():
        r, c = (differences[stratum, key] for key in ("R_at_P0_C0", "C_at_P0_R0"))
        assert f"| {label} | {r['families']} | {r['F1_pp']:+.3f} | {c['F1_pp']:+.3f} |" in text
