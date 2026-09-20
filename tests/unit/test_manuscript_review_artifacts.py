from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.parse import unquote, urlsplit

from benchmark_tools.prepare_ob_candidate_neighborhood import record

ROOT = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


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


def test_all_review_figures_link_to_existing_full_resolution_sources():
    parsed = Images()
    parsed.feed((ROOT / "PUBLICATION_MANUSCRIPT_REVIEW_20260920.html").read_text())
    assert len(parsed.images) == 8
    for source, link, alt in parsed.images:
        assert source == link and alt
        assert (ROOT / unquote(urlsplit(source).path)).is_file()


def test_review_output_and_figure_identities():
    report = json.loads((ROOT / "manuscript_local_assets_review_20260920.json").read_text())
    image_urls = {row["url"] for row in report["occurrences"] if row["kind"] == "Image"}
    figures = [item for item in report["targets"]
               if item["path"].split("/benchmark_tools/results/", 1)[-1] in image_urls]
    assert len(figures) == 8
    # The dated review remains fixed; the living manuscript and notes can evolve.
    for item in [report["html"], *figures]:
        relative = item["path"].split("/benchmark_tools/results/", 1)[-1]
        current = record(ROOT / relative)
        assert (current["sha256"], current["bytes"]) == (item["sha256"], item["bytes"])
    assert report["unique_targets"] == 162
    assert report["local_occurrences"] == 180
    assert report["publication_ready"] is False
