from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.render_publication_bibliography import audit, entries

ROOT = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def test_v4_changes_only_reviewed_matplotlib_journal_name():
    old = json.loads((ROOT / "publication_bibliography_20260919_v3.csl.json").read_text())
    new = json.loads((ROOT / "publication_bibliography_20260920_v4.csl.json").read_text())
    expected = deepcopy(old)
    row, = [r for r in expected if r["id"] == "matplotlib2007"]
    assert row["DOI"] == "10.1109/MCSE.2007.55"
    assert row["container-title"] == "Computing in Science &amp; Engineering"
    row["container-title"] = "Computing in Science & Engineering"
    assert len(new) == 37 and new == expected


def test_current_render_contains_corrected_text_and_all_fields():
    records = json.loads((ROOT / "publication_bibliography_20260920_v4.csl.json").read_text())
    document = json.loads((ROOT / "publication_bibliography_review_20260920/bibliography.pandoc.json").read_text())
    assert audit(records, document)["all_checked_fields_present"] is True
    text = entries(document)["matplotlib2007"]
    assert "Computing in Science & Engineering" in text
    assert "&amp;" not in text


@pytest.mark.parametrize("name,field,folder", [
    ("publication_bibliography_visual_review_20260920.json", "review_artifacts", "publication_bibliography_review_20260920"),
    ("publication_citation_fields_review_20260920.json", "artifacts", "publication_bibliography_review_20260920_v5")])
def test_visual_review_records_bind_committed_artifacts(name, field, folder):
    report = json.loads((ROOT / name).read_text())
    for item in [report["input"], report["output"], *report[field]]:
        path = ROOT / Path(item["path"]).name
        if item in report[field]:
            path = ROOT / folder / path.name
        observed = record(path)
        assert (observed["sha256"], observed["bytes"]) == (item["sha256"], item["bytes"])


def test_v5_changes_only_two_publisher_reviewed_fields():
    old = json.loads((ROOT / "publication_bibliography_20260920_v4.csl.json").read_text())
    new = json.loads((ROOT / "publication_bibliography_20260920_v5.csl.json").read_text())
    expected = deepcopy(old)
    indexed = {r["id"]: r for r in expected}
    assert indexed["fastme2015"]["title"].endswith(": Table 1.")
    assert indexed["Li2006TreeFam"]["issue"] == "90001"
    indexed["fastme2015"]["title"] = "FastME 2.0: A Comprehensive, Accurate, and Fast Distance-Based Phylogeny Inference Program"
    indexed["Li2006TreeFam"]["issue"] = "suppl_1"
    assert len(new) == 37 and new == expected


def test_v5_render_preserves_inventory_and_corrected_issue():
    records = json.loads((ROOT / "publication_bibliography_20260920_v5.csl.json").read_text())
    document = json.loads((ROOT / "publication_bibliography_review_20260920_v5/bibliography.pandoc.json").read_text())
    assert audit(records, document)["all_checked_fields_present"] is True
    rendered = entries(document)
    assert "Table 1" not in rendered["fastme2015"]
    assert "Issue: suppl_1" in rendered["Li2006TreeFam"]
