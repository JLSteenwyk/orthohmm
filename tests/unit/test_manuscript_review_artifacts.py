from html.parser import HTMLParser
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote, urlsplit

import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record

ROOT = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def test_null_score_main_review_retains_closed_artifacts_and_scope():
    review = json.loads((ROOT / "publication_main_null_visual_review_20261002.json").read_bytes())
    repo = ROOT.parents[1]
    assert review["render_time_source_commit"] == "6ac1ccbb4bd4fc793333ceddb700329124757a76"
    assert review["page_count"] == 8 and review["all_eight_pages_inspected"] is True
    assert review["bounds_violations"] == 0 and review["observed_clipping_or_overlap"] is False
    assert review["null_results_page"] == 4 and review["null_protocol_page"] == 2
    assert review["local_occurrences"] == 41 and review["unique_targets"] == 40
    assert len(review["citation_ids"]) == 16 and not review["untracked_targets"]
    assert len(review["closed_review_artifacts"]) == 13
    for item in review["closed_review_artifacts"]:
        actual = record(repo / item["path"])
        assert (actual["bytes"], actual["sha256"]) == (item["bytes"], item["sha256"])
    for key in ("native_inference_or_scoring_rerun", "benchmark_scores_or_defaults_changed",
                "controlled_timing_executed", "new_study_archive_built", "data_rights_cleared",
                "public_release_or_deposition_executed", "publication_ready", "other_linked_figures_newly_revalidated"):
        assert review[key] is False
    assert review["new_null_figure_separately_inspected"] is True
    html = (ROOT / "publication_main_review_20261002_v2.html").read_text()
    parts = []
    class Text(HTMLParser):
        def handle_data(self, data):
            parts.append(data)
    Text().feed(html)
    visible = " ".join(" ".join(parts).split())
    assert "89.83%, 100% and 100%" in visible
    assert "not a predicted ortholog or an observed pipeline false positive" in visible
    assert "90,000 independent pairs, not 180,000 independent observations" in visible
    assert "No coefficients, thresholds or defaults were fitted or promoted" in visible
    assert "No prefilter or biological inference ran" in visible
    assert "Not submission-ready" in visible
    assert "figures_frozen_null_scores_20261002_v2/frozen_null_scores.pdf" in html


def test_october_2_main_review_retains_closed_artifact_identities():
    review = json.loads((ROOT / "publication_main_visual_review_20261002.json").read_text())
    repo = ROOT.parents[1]
    items = [review[k] for k in ("render", "html", "print", "pdf", "bounds_review")]
    items.extend(review["reviewed_pages"])
    assert len(items) == 12 and len(review["reviewed_pages"]) == 7
    for item in items:
        actual = record(repo / item["path"])
        assert (actual["bytes"], actual["sha256"]) == (item["bytes"], item["sha256"])
    assert review["render_time_source_commit"] == "cb5bda8083c24d4915bab5dd78de6d9b982195be"
    assert review["page_count"] == 7 and review["bounds_violations"] == 0
    assert review["local_link_occurrences"] == 37 and review["unique_local_targets"] == 36
    assert len(review["citation_ids"]) == 16 and review["untracked_targets"] == []
    assert review["all_seven_pages_inspected"] is True
    assert review["linked_figures_newly_visually_revalidated"] is False
    assert review["new_review_archive_built"] is review["publication_ready"] is False
    html = (ROOT / "publication_main_review_20261002.html").read_text()
    parts = []
    class Text(HTMLParser):
        def handle_data(self, data):
            parts.append(data)
    Text().feed(html)
    visible = " ".join(" ".join(parts).split())
    assert "updated 2 October 2026" in visible and "Not submission-ready" in visible
    for link in ("SWISS_DESCRIPTIVE_COMPONENT_20261002.md", "SWISS_RAW_ARCHIVE_RESTORATION_20261002.md"):
        assert link in html and (ROOT / link).is_file()
    # Living sources can evolve; this assertion binds only the dated closed render.
    assert "private archives remain unuploaded" in visible


def test_v10_opening_retains_linked_bibliography_history():
    html = ROOT / "PUBLICATION_MANUSCRIPT_REVIEW_20260925_v10.html"
    text = html.read_text()
    assert "25 September 2026" in text
    assert "PUBLICATION_MANUSCRIPT_BIBLIOGRAPHY_PROVENANCE_20260925.md" in text
    assert "FastME table-label title suffix" not in text
    note = (ROOT / "PUBLICATION_MANUSCRIPT_BIBLIOGRAPHY_PROVENANCE_20260925.md").read_text()
    assert "FastME table-label title suffix" in note
    report = json.loads((ROOT / "manuscript_asset_review_20260925_v10.json").read_text())
    assert report["html"]["sha256"] == record(html)["sha256"]
    assert report["unique_targets"] == 175 and report["local_occurrences"] == 189
    assert report["untracked_targets"] == []
    assert hashlib.sha256((ROOT / "PUBLICATION_MANUSCRIPT_REVIEW_20260925_v10.pdf").read_bytes()).hexdigest() == "ee70c37471bc07d5b78fc022d229e4e6e7e9741e89f65442b61b646e71e39bbc"


def test_v10_retains_three_caption_groups():
    parsed = FigureGroups()
    parsed.feed((ROOT / "PUBLICATION_MANUSCRIPT_REVIEW_20260925_v10.html").read_text())
    assert len(parsed.groups) == 3
    for group in parsed.groups:
        assert len(group["images"]) == 1
        assert "break-inside: avoid" in group["style"]
        assert "Supplementary Figure:" in "".join(group["text"])


class FigureGroups(HTMLParser):
    def __init__(self):
        super().__init__()
        self.groups = []
        self.active = None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "div" and "review-figure" in attrs.get("class", "").split():
            assert self.active is None
            self.active = dict(style=attrs["style"], images=[], text=[])
        if tag == "img" and self.active is not None:
            self.active["images"].append(attrs["src"])

    def handle_data(self, text):
        if self.active is not None:
            self.active["text"].append(text)

    def handle_endtag(self, tag):
        if tag == "div" and self.active is not None:
            self.groups.append(self.active)
            self.active = None


def test_v9_supplementary_figures_keep_explicit_captions():
    parsed = FigureGroups()
    parsed.feed((ROOT / "PUBLICATION_MANUSCRIPT_REVIEW_20260923_v9.html").read_text())
    assert len(parsed.groups) == 3
    for group, title in zip(parsed.groups, ("Duplication Annotations", "Sequence Identity", "Fragment Annotations")):
        assert len(group["images"]) == 1
        assert "break-inside: avoid" in group["style"]
        text = " ".join("".join(group["text"]).split())
        assert "Supplementary Figure: " + title + "." in text
    report = json.loads((ROOT / "manuscript_asset_review_20260923_v9.json").read_text())
    actual = record(ROOT / "PUBLICATION_MANUSCRIPT_REVIEW_20260923_v9.html")
    assert actual["sha256"] == report["html"]["sha256"]
    assert report["unique_targets"] == 180 and report["local_occurrences"] == 198


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
