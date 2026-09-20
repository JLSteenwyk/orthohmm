from copy import deepcopy
import json
from pathlib import Path

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


def test_visual_review_records_bind_committed_artifacts():
    report = json.loads((ROOT / "publication_bibliography_visual_review_20260920.json").read_text())
    for item in [report["input"], report["output"], *report["review_artifacts"]]:
        path = ROOT / Path(item["path"]).name
        if item in report["review_artifacts"]:
            path = ROOT / "publication_bibliography_review_20260920" / path.name
        observed = record(path)
        assert (observed["sha256"], observed["bytes"]) == (item["sha256"], item["bytes"])
