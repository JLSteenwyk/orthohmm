import json

import fitz
import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.review_manuscript_pdf import review


def fixture(tmp_path):
    pdf = tmp_path / "document.pdf"
    with fitz.open() as document:
        for i in range(3):
            page = document.new_page()
            page.insert_text((50, 50), "selected paragraph" if i == 1 else "other page")
        document.save(pdf)
    html = tmp_path / "source.html"
    html.write_text("<p>selected paragraph</p>")
    assets = tmp_path / "assets.json"
    assets.write_text(json.dumps(dict(status="manuscript_review_rendered", html=record(html),
        sources=[], targets=[], local_occurrences=0, unique_targets=0, untracked_targets=[])))
    return pdf, assets, html


def test_selects_matching_page_and_neighbors(tmp_path):
    pdf, assets, _ = fixture(tmp_path)
    result = review(pdf, assets, tmp_path / "output", ["selected paragraph"])
    assert result["phrase_matches"] == {"selected paragraph": [2]}
    assert len(result["rendered_pages"]) == 3
    assert result["bounds_violations"] == []
    assert result["visual_review_complete"] is False
    with pytest.raises(FileExistsError):
        review(pdf, assets, tmp_path / "output", ["selected paragraph"])


def test_reject_missing_selector(tmp_path):
    pdf, assets, _ = fixture(tmp_path)
    with pytest.raises(ValueError, match="absent"):
        review(pdf, assets, tmp_path / "output", ["not present"])
    assert not (tmp_path / "output").exists()


def test_reject_changed_source(tmp_path):
    pdf, assets, html = fixture(tmp_path)
    html.write_text("changed")
    with pytest.raises(ValueError):
        review(pdf, assets, tmp_path / "output", ["selected paragraph"])
