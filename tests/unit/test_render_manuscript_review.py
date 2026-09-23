import json
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
    with pytest.raises(FileExistsError):
        render(tmp_path, draft, output, report)


def test_requires_sibling_html_and_distinct_outputs(tmp_path):
    draft = tmp_path / "draft.md"
    for output, report in [(tmp_path / "elsewhere/review.html", tmp_path / "report.json"),
                           (tmp_path / "same", tmp_path / "same")]:
        with pytest.raises(ValueError, match="distinct review"):
            render(tmp_path, draft, output, report)
