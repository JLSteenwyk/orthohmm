from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.parse import unquote, urlsplit

import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record

ROOT = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def test_v6_updated_results_and_pdf_identity():
    html = ROOT / "PUBLICATION_MANUSCRIPT_REVIEW_20260923_v6.html"
    text = html.read_text()
    assert "23 September 2026" in text
    assert "276 genes" in text
    assert "found14SwissTrees" not in text
    assert "QFO_VATB_PARTITION_TRACE_20260923.md" in text
    report = json.loads((ROOT / "manuscript_asset_review_20260923_v6.json").read_text())
    assert report["unique_targets"] == 178
    assert report["local_occurrences"] == 196
    assert report["untracked_targets"] == []
    review = json.loads((ROOT / "manuscript_pdf_review_20260923_v6.json").read_text())
    assert review["embedded_images"] == 11
    assert review["visually_reviewed_pages"] == [12, 17]
    assert review["publication_ready"] is False
    for item in (report["html"], review["pdf"], review["html"], review["asset_report"]):
        actual = record(ROOT / Path(item["path"]).name)
        assert (actual["bytes"], actual["sha256"]) == (item["bytes"], item["sha256"])


class Images(HTMLParser):
    def __init__(self):
        super().__init__()
        self.link = None
        self.images = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "a":
            self.link = attrs.get("href")
        elif tag == "img":
            self.images.append((attrs.get("src"), self.link, attrs.get("alt")))

    def handle_endtag(self, tag):
        if tag == "a":
            self.link = None


@pytest.mark.parametrize("date", ["20260920", "20260923"])
def test_all_review_figures_link_to_existing_full_resolution_sources(date):
    parsed = Images()
    parsed.feed((ROOT / f"PUBLICATION_MANUSCRIPT_REVIEW_{date}.html").read_text())
    assert len(parsed.images) == 8
    for source, link, alt in parsed.images:
        assert source == link and alt
        assert (ROOT / unquote(urlsplit(source).path)).is_file()


@pytest.mark.parametrize("date,targets,occurrences", [("20260920", 162, 180), ("20260923", 165, 183)])
def test_review_output_and_figure_identities(date, targets, occurrences):
    report = json.loads((ROOT / f"manuscript_local_assets_review_{date}.json").read_text())
    image_urls = {row["url"] for row in report["occurrences"] if row["kind"] == "Image"}
    figures = [item for item in report["targets"]
               if item["path"].split("/benchmark_tools/results/", 1)[-1] in image_urls]
    assert len(figures) == 8
    # The dated review remains fixed; the living manuscript and notes can evolve.
    for item in [report["html"], *figures]:
        relative = item["path"].split("/benchmark_tools/results/", 1)[-1]
        current = record(ROOT / relative)
        assert (current["sha256"], current["bytes"]) == (item["sha256"], item["bytes"])
    assert report["unique_targets"] == targets
    assert report["local_occurrences"] == occurrences
    assert report["publication_ready"] is False
