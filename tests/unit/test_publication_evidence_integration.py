"""Check revised manuscript claims against retained evidence, not new scoring."""

import hashlib
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"
OLD = (RESULTS / "PUBLICATION_MAIN_TEXT_20261004_v2.md").read_text()
NEW = (RESULTS / "PUBLICATION_MAIN_TEXT_20261004_v3.md").read_text()


def read(path):
    return json.loads((RESULTS / path).read_bytes())


@pytest.mark.parametrize("path,expected", [
    ("PUBLICATION_MAIN_TEXT_20261004_v2.md", "8c4635d8bb5041937ffc60f3e5155d6e60e6a53b808d174c78d6801e2592d337"),
    ("publication_main_with_figures_20261004_v2/document.pdf", "07f2411b3e59f8295c928b8164b30ee6ba0c76cc4a470121871941200a73d1f2"),
])
def test_previous_source_and_render_are_unchanged(path, expected):
    assert hashlib.sha256((RESULTS / path).read_bytes()).hexdigest() == expected


def test_abstract_results_and_references_preserve_original_scientific_text():
    for heading, end in [("## Abstract\n", "## Methods\n"), ("## References\n", None)]:
        original = OLD.split(heading, 1)[1]
        current = NEW.split(heading, 1)[1]
        if end:
            original, current = original.split(end, 1)[0], current.split(end, 1)[0]
        assert original == current
    old_results = OLD.split("## Results\n", 1)[1].split("## Discussion And Limitations\n", 1)[0]
    new_results = NEW.split("## Results\n", 1)[1].split("## Discussion And Limitations\n", 1)[0]
    assert new_results.startswith(old_results)


def test_inventory_claims_match_actual_families_and_partition_roles():
    evidence = read("development_family_inventory_20261004/inventory.json")
    assert evidence["json_files_scanned"] == 1899
    assert sum(f["origin"] == "git_snapshot" for f in evidence["files"]) == 1785
    assert sum(f["origin"] == "frozen_local_only" for f in evidence["files"]) == 114
    assert evidence["summary"]["OrthoBench"]["scored_blocks"] == 134
    assert evidence["summary"]["OrthoBench"]["scored_files"] == 84
    assert evidence["summary"]["OrthoBench"]["family_block_associations"] == 7770
    assert evidence["summary"]["QfO_SwissTrees"]["scored_blocks"] == 77
    assert evidence["summary"]["QfO_SwissTrees"]["scored_files"] == 16
    assert evidence["summary"]["QfO_SwissTrees"]["family_block_associations"] == 1386
    for row in evidence["family_rows"]:
        assert row["scored_evidence_blocks"] > 0
        if row["dataset"] != "OrthoBench":
            continue
        role = row["original_partition"]
        assert row["reported_partition_block_counts"][role] == (28 if role == "development" else 18)
        assert row["reported_partition_block_counts"]["all"] == 21
    for phrase in ["18 declared validation blocks", "28 development blocks", "21 all-partition blocks",
                   "134 OrthoBench score blocks", "77 SwissTrees blocks", "7,770 and 1,386"]:
        assert phrase in NEW
    assert "not untouched publication\nvalidation" in NEW
    assert "not a\ncomplete causal tuning history" in NEW


@pytest.mark.parametrize("case,identity", [("baseline", 8), ("cluster_order_only", 0), ("one_ulp_scores_only", 7)])
def test_unattached_satellite_values_are_identifiers_not_counts(case, identity):
    fixture = read("candidate_trace_variation_20261004/diagnostic.json")["fixture"]
    row = next(c for c in fixture["cases"] if c["case"] == case)
    assert row["unattached_satellites"] == [identity]
    assert row["merges"] == 8 and row["iterations"] == 2
    assert sum(len(group) == 1 for group in row["partition"]) == 1
    assert "satellite IDs 8, 0 and 7" in NEW
    assert "one unattached\nsatellite and eight accepted merges in every case" in NEW


def test_trace_counts_and_support_difference_match_retained_evidence():
    report = read("candidate_trace_variation_20261004/diagnostic.json")
    counts = [point["comparison"]["common_semantic_merges"] for point in report["points"]]
    assert counts == [8428, 8439, 8426]
    assert "8,428, 8,439 and 8,426" in NEW
    for point in report["points"]:
        comparison = point["comparison"]
        assert comparison["original_merges"] == comparison["native_merges"] == 8440
        assert comparison["common_source_cluster_ids_changed"] == comparison["common_target_cluster_ids_changed"] == 0
        assert comparison["common_feature_differences"]["support"]["maximum_absolute_delta"] == 2.1316282072803006e-14
    assert "2.1316282072803006e-14" in NEW
    assert "does not establish historical relabeling, score-bit provenance" in NEW


def test_native_costs_are_descriptive_and_do_not_replace_original_cached_costs():
    report = read("factorial_native_resource_linkage_20261004/linkage.json")
    for row in report["summaries"]:
        minutes = row["resources"]["wall_seconds"]["median"] / 60
        assert f"{minutes:.3f}" in NEW
    assert [row["exact_partition_repeats"] for row in report["summaries"]] == [3, 1]
    assert "patched native runtime differs from that deployment" in NEW
    assert "costs of the original\ncached executions" in NEW
    assert "Affected-group gene counts are not counts of genes moved" in NEW


def test_selected_qfo_stage_associations_preserve_unavailable_full_costs():
    report = read("qfo_orthohmm_stage_metadata_20261004/register.json")
    assert report["method_dataset_cells"] == 24
    selected = [row["qfo_stage_provenance"] for row in report["rows"] if "qfo_stage_provenance" in row]
    assert {row["cell"] for row in selected} == {"p1_c0_r0", "p1_c1_r1"}
    for row in selected:
        assert row["full_pipeline_wall_s"] is None
        assert row["full_pipeline_cpu_s"] is None
        assert row["full_pipeline_peak_memory_bytes"] is None
        assert row["observations_are_not_summed"] is True
    assert "one observation, not two repeats" in NEW
    assert "complete cached-execution\ncosts remain unknown" in NEW


@pytest.mark.parametrize("name", [
    "DEVELOPMENT_FAMILY_INVENTORY_RESULT_20261004.md",
    "DEVELOPMENT_FAMILY_INVENTORY_PROTOCOL_20261004.md",
    "CANDIDATE_TRACE_VARIATION_RESULT_20261004.md",
    "FACTORIAL_NATIVE_RESOURCE_LINKAGE_RESULT_20261004.md",
    "QFO_ORTHOHMM_STAGE_LINKAGE_RESULT_20261004.md",
    "FACTORIAL_RESOURCE_RESULT_20261004.md",
    "ALL_BENCHMARK_METADATA_INTEGRATION_RESULT_20261004.md",
])
def test_new_evidence_links_exist(name):
    assert f"]({name})" in NEW and (RESULTS / name).is_file()


def test_new_source_is_not_claimed_as_rendered_or_archived():
    assert "Neither is a rendered or archived copy of this revised Markdown" in NEW
    assert "those steps have not been executed here" in NEW
    assert "No new endpoint selection, score, method default or independent-validation" in NEW
