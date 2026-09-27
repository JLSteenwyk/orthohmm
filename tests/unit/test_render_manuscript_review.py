import json
from html.parser import HTMLParser
from pathlib import Path
import shutil
import subprocess

import pytest

from benchmark_tools.render_manuscript_review import local_assets, render


def node(kind, url):
    return {"t": kind, "c": [["", [], []], [], [url, ""]]}


def test_local_inventory_decodes_deduplicates_and_preserves_occurrences(tmp_path):
    asset = tmp_path / "figure name.png"
    asset.write_bytes(b"placeholder")
    document = [node("Link", "figure%20name.png#anchor"), node("Image", "figure%20name.png"),
                node("Link", "#section"), node("Link", "https://example.org/"),
                node("Link", "//example.org/asset")]
    rows, targets = local_assets(document, tmp_path / "draft.md", tmp_path)
    assert len(rows) == 2 and len(targets) == 1
    assert [r["kind"] for r in rows] == ["Link", "Image"]
    assert targets["figure name.png"]["bytes"] == 11


@pytest.mark.parametrize("url", ["../escape", "/absolute", "%2e%2e/escape"])
def test_escape_rejected(tmp_path, url):
    with pytest.raises(ValueError, match="escapes"):
        local_assets([node("Link", url)], tmp_path / "draft.md", tmp_path)


def test_missing_asset_rejected(tmp_path):
    with pytest.raises(FileNotFoundError):
        local_assets([node("Image", "missing.png")], tmp_path / "draft.md", tmp_path)


def test_symlink_escape_rejected(tmp_path):
    (tmp_path / "link").symlink_to(tmp_path.parent)
    with pytest.raises(ValueError, match="escapes"):
        local_assets([node("Link", "link/file")], tmp_path / "draft.md", tmp_path)


@pytest.mark.skipif(shutil.which("pandoc") is None, reason="pandoc required")
def test_real_pandoc_render_and_no_overwrite(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    draft = tmp_path / "draft.md"
    draft.write_text("# Review\n\n[Evidence](evidence.json)\n")
    (tmp_path / "evidence.json").write_text("{}\n")
    subprocess.run(["git", "add", "draft.md", "evidence.json"], cwd=tmp_path, check=True)
    output, report = tmp_path / "review.html", tmp_path / "review.json"
    result = render(tmp_path, draft, output, report)
    assert result == json.loads(report.read_text())
    assert result["local_occurrences"] == result["unique_targets"] == 1
    assert result["untracked_targets"] == []
    assert result["publication_ready"] is False
    assert 'href="evidence.json"' in output.read_text()
    class Headings(HTMLParser):
        count = 0

        def handle_starttag(self, tag, attrs):
            if tag == "h1":
                self.count += 1

    headings = Headings()
    headings.feed(output.read_text())
    assert headings.count == 1
    assert "<title>OrthoHMM publication working draft</title>" in output.read_text()
    assert "@page { margin: 18mm; }" in output.read_text()
    assert any(Path(s["path"]).name == "manuscript_review_print.html" for s in result["sources"])
    with pytest.raises(FileExistsError):
        render(tmp_path, draft, output, report)


def test_requires_sibling_html_and_distinct_outputs(tmp_path):
    draft = tmp_path / "draft.md"
    for output, report in [(tmp_path / "elsewhere/review.html", tmp_path / "report.json"),
                           (tmp_path / "same", tmp_path / "same")]:
        with pytest.raises(ValueError, match="distinct review"):
            render(tmp_path, draft, output, report)


@pytest.mark.skipif(shutil.which("pandoc") is None, reason="pandoc required")
@pytest.mark.parametrize("case", ["valid", "missing", "duplicate", "no_bibliography", "escape"])
def test_citation_rendering_and_rejection(tmp_path, case):
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    draft = repo / "draft.md"
    draft.write_text("# Review\n\nMethod [@method].\n\n## References\n")
    bibliography = repo / "references.json"
    entry = {"id": "method", "type": "article-journal", "title": "Test Method",
             "author": [{"family": "Example"}], "issued": {"date-parts": [[2020]]}}
    entries = [entry]
    if case == "duplicate":
        entries.append(entry)
    elif case == "missing":
        entries = []
    elif case == "escape":
        bibliography = tmp_path / "references.json"
    bibliography.write_text(json.dumps(entries))
    output, report = repo / "review.html", repo / "review.json"
    if case != "valid":
        with pytest.raises(ValueError):
            render(repo, draft, output, report, None if case == "no_bibliography" else bibliography)
        assert not output.exists() and not report.exists()
        return
    result = render(repo, draft, output, report, bibliography)
    assert result["citation_ids"] == ["method"]
    assert result["stderr"] == {"parse": "", "render": ""}
    assert 'id="ref-method"' in output.read_text()
    assert "Test Method" in output.read_text()
    assert any(Path(s["path"]) == bibliography for s in result["sources"])
