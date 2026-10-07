"""The new review removes a mutable asset without changing scientific content."""

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "benchmark_tools/results"


def test_review_asset_revision_is_exactly_the_two_declared_edits():
    receipt = json.loads((BASE / "publication_main_review_asset_revision_20261007_v3.json").read_text())
    for key in ("parent", "output"):
        ref = receipt[key]
        data = (ROOT / ref["path"]).read_bytes()
        assert len(data) == ref["bytes"]
        assert hashlib.sha256(data).hexdigest() == ref["sha256"]
    parent = (ROOT / receipt["parent"]["path"]).read_text()
    expected = parent
    assert len(receipt["edits"]) == 2
    for edit in receipt["edits"]:
        assert expected.count(edit["from"]) == edit["occurrences"] == 1
        expected = expected.replace(edit["from"], edit["to"], 1)
    actual = (ROOT / receipt["output"]["path"]).read_text()
    assert actual == expected
    assert "[progress ledger](PUBLICATION_PROGRESS.md)" not in actual
    assert "progress ledger records completed work and unmet" in actual
    assert [l for l in actual.splitlines() if l.startswith("|")] == [
        l for l in parent.splitlines() if l.startswith("|")]
    for key in ("old_render_or_pdf_modified", "benchmark_or_scoring_reexecuted",
                "new_render_or_package_proved", "publication_ready"):
        assert receipt[key] is False


def test_failed_v2_asset_review_is_retained_not_relabelled_successful():
    revision = json.loads((BASE / "publication_main_review_asset_revision_20261007_v3.json").read_text())
    path = ROOT / revision["failed_review"]["path"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == revision["failed_review"]["sha256"]
    failure = json.loads(path.read_text())
    assert failure["returncode"] == 1 and failure["error_type"] == "ValueError"
    assert failure["status"] == "pdf_asset_review_failed_before_rasterization"
    assert "PUBLICATION_PROGRESS.md" in failure["error"]
    assert failure["visual_review_complete"] is False and failure["publication_ready"] is False
    assert not (ROOT / failure["requested_output"]).exists()
    render = json.loads((ROOT / failure["assets"]).read_text())
    binding = next(r for r in render["targets"] if r["path"].endswith("PUBLICATION_PROGRESS.md"))
    assert binding["sha256"] == failure["recorded_ledger_sha256"]
    assert binding["bytes"] == failure["recorded_ledger_bytes"]
    print_receipt = json.loads((BASE / "publication_main_print_20261007_v2/print.json").read_text())
    pdf = ROOT / failure["pdf"]
    assert hashlib.sha256(pdf.read_bytes()).hexdigest() == print_receipt["pdf"]["sha256"]
    assert print_receipt["page_count"] == 21
