"""Pinned presentation replay, relocation and actual final-inventory checks."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import fitz
import pytest

from benchmark_tools.assemble_publication_review import assemble


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"


def ref(root, path):
    content = path.read_bytes()
    return {"path": str(path.relative_to(root)), "bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest()}


def write_json(path, value):
    path.write_text(json.dumps(value))
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def selection(tmp_path):
    root = tmp_path / "inputs"
    root.mkdir()
    with fitz.open() as main:
        for index in range(2):
            page = main.new_page()
            page.insert_text((40, 40), "Main page %d" % (index + 1))
        main[0].insert_link({"kind": fitz.LINK_URI, "from": fitz.Rect(40, 50, 100, 70),
                             "uri": "file:///historical/figure.pdf"})
        main[1].insert_link({"kind": fitz.LINK_URI, "from": fitz.Rect(40, 50, 100, 70),
                             "uri": "file:///historical/nonfigure.json"})
        main[1].insert_link({"kind": fitz.LINK_URI, "from": fitz.Rect(40, 80, 100, 100),
                             "uri": "https://example.org/reference"})
        main.save(root / "main.pdf")
    figures = []
    for number in (1, 2):
        path = root / ("figure_%d.pdf" % number)
        with fitz.open() as figure:
            page = figure.new_page(width=320 + number, height=210 + number)
            page.insert_text((40, 40), "Figure %d" % number)
            figure.save(path)
        figures.append({"number": number, "title": "Figure %d" % number,
                        "caption": "Retained descriptive figure; not inference.",
                        "pdf": ref(root, path),
                        "link_aliases": ["/historical/figure.pdf"] if number == 2 else []})
    chosen = {"schema": "publication_figure_selection_v1", "main": ref(root, root / "main.pdf"),
              "main_pages": 2, "figures": figures, "provenance": []}
    path = root / "selection.json"
    anchor = write_json(path, chosen)
    return root, path, anchor, chosen


@pytest.mark.parametrize("optimized", [False, True])
def test_relocated_cli_preserves_sources_and_actions(selection, tmp_path, optimized):
    root, selection_path, anchor, chosen = selection
    runner = root / "assemble.py"
    shutil.copyfile(ROOT / "benchmark_tools/assemble_publication_review.py", runner)
    command = [sys.executable, "-I", "-B", *(["-O"] if optimized else []), str(runner),
               "--root", str(root), "--selection", str(selection_path), "--selection-sha256", anchor,
               "--output", str(tmp_path / "assembled")]
    completed = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert completed.stderr == ""
    report = json.loads((tmp_path / "assembled/assembly.json").read_text())
    assert (report["main_pages"], report["guide_pages"], report["figure_pages"], report["total_pages"]) == (2, 1, 2, 5)
    assert report["publication_ready"] is report["visual_review_complete"] is False
    assert report["redirected_main_figure_links"] == [{"main_page": 1,
        "source": "/historical/figure.pdf", "combined_page": 5}]
    assert report["restored_original_file_uri_actions"] == [{"main_page": 2,
        "uri": "file:///historical/nonfigure.json"}]
    with fitz.open(tmp_path / "assembled/document.pdf") as combined:
        assert combined.get_toc() == [[1, "Main Text", 1], [1, "Figure Guide", 3],
                                     [1, "A1. Figure 1", 4], [1, "A2. Figure 2", 5]]
        assert combined[0].get_links()[0]["page"] == 4
        assert combined[1].get_links()[1]["uri"] == "https://example.org/reference"
        assert [link["page"] for link in combined[2].get_links()] == [3, 4]
        for item in report["preserved_source_pages"]:
            with fitz.open(item["source"]["path"]) as source:
                before, after = source[item["original_page"] - 1], combined[item["combined_page"] - 1]
                assert before.rect == after.rect
                assert before.get_text("words") == after.get_text("words")
                assert before.get_pixmap().samples == after.get_pixmap().samples


@pytest.mark.parametrize("damage", ["selection", "pdf", "unsafe", "duplicate", "alias", "page_count", "caption"])
def test_bad_inputs_fail_before_output(selection, tmp_path, damage):
    root, path, anchor, chosen = selection
    if damage == "selection":
        path.write_text(path.read_text() + " ")
    elif damage == "pdf":
        target = root / chosen["figures"][0]["pdf"]["path"]
        target.write_bytes(target.read_bytes() + b"corrupt")
    else:
        if damage == "unsafe":
            chosen["main"]["path"] = "../main.pdf"
        elif damage == "duplicate":
            chosen["figures"][1]["pdf"] = chosen["figures"][0]["pdf"]
        elif damage == "alias":
            chosen["figures"][0]["link_aliases"] = ["/historical/figure.pdf"]
        elif damage == "page_count":
            chosen["main_pages"] = 9
        elif damage == "caption":
            chosen["figures"][0]["caption"] = "Oversized caption " * 1000
        anchor = write_json(path, chosen)
    with pytest.raises(ValueError):
        assemble(root, path, anchor, tmp_path / "bad")
    assert not (tmp_path / "bad").exists()


def test_existing_output_untouched(selection, tmp_path):
    root, path, anchor, _ = selection
    output = tmp_path / "existing"
    output.mkdir()
    marker = output / "marker"
    marker.write_text("preserve")
    with pytest.raises(FileExistsError, match="Refusing existing"):
        assemble(root, path, anchor, output)
    assert list(output.iterdir()) == [marker]
    assert marker.read_text() == "preserve"


def test_symlink_escape_rejected(selection, tmp_path):
    root, path, _, chosen = selection
    outside = tmp_path / "outside.pdf"
    shutil.copyfile(root / "main.pdf", outside)
    (root / "escape.pdf").symlink_to(outside)
    chosen["main"]["path"] = "escape.pdf"
    anchor = write_json(path, chosen)
    with pytest.raises(ValueError, match="escapes root"):
        assemble(root, path, anchor, tmp_path / "bad")
    assert not (tmp_path / "bad").exists()


def test_actual_final_inventory_and_closed_artifacts():
    chosen = json.loads((RESULTS / "publication_figure_selection_20261004.json").read_text())
    report = json.loads((RESULTS / "publication_main_with_figures_20261004/assembly.json").read_text())
    assert report["main_pages"] == 12
    assert report["guide_pages"] == 4
    assert report["figure_pages"] == 16
    assert report["total_pages"] == 32
    assert len(report["preserved_source_pages"]) == 28
    assert len(report["rendered_pages"]) == 20
    assert chosen["figures"][14]["title"] == "All-method OrthoBench error strata"
    assert chosen["figures"][15]["title"] == "Completed shared-host resource panel"
    assert "24 eligible" in chosen["figures"][15]["caption"]
    assert "Three cells" in chosen["figures"][15]["caption"]
    assert "potentially tool-dependent contention" in chosen["figures"][15]["caption"]
    assert report["publication_ready"] is False
    for reference in [report["pdf"], *report["inputs"], *report["rendered_pages"]]:
        content = Path(reference["path"]).read_bytes()
        assert (len(content), hashlib.sha256(content).hexdigest()) == (reference["bytes"], reference["sha256"])
    with fitz.open(report["pdf"]["path"]) as combined:
        assert len(combined.get_toc()) == 18
        assert combined[31].get_text("words")
        assert sum(len(combined[index].get_links()) for index in range(12, 16)) == 16
        for item in report["preserved_source_pages"]:
            with fitz.open(item["source"]["path"]) as source:
                before, after = source[item["original_page"] - 1], combined[item["combined_page"] - 1]
                assert before.rect == after.rect
                assert before.get_text("words") == after.get_text("words")
                assert before.get_pixmap().samples == after.get_pixmap().samples


def test_actual_inventory_replays_without_checkout(tmp_path):
    selection_path = RESULTS / "publication_figure_selection_20261004.json"
    chosen = json.loads(selection_path.read_text())
    root = tmp_path / "copied"
    references = [chosen["main"], *chosen["provenance"],
                  *[item["pdf"] for item in chosen["figures"]]]
    for reference in references:
        destination = root / reference["path"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / reference["path"], destination)
    copied_selection = root / "selection.json"
    shutil.copyfile(selection_path, copied_selection)
    runner = root / "assemble.py"
    shutil.copyfile(ROOT / "benchmark_tools/assemble_publication_review.py", runner)
    anchor = hashlib.sha256(copied_selection.read_bytes()).hexdigest()
    completed = subprocess.run([sys.executable, "-I", "-B", "-O", str(runner), "--root", str(root),
        "--selection", str(copied_selection), "--selection-sha256", anchor,
        "--output", str(tmp_path / "replayed")], cwd=tmp_path, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr
    assert completed.stderr == ""
    replayed = json.loads((tmp_path / "replayed/assembly.json").read_text())
    assert all(Path(reference["path"]).is_relative_to(root) for reference in replayed["inputs"])
    assert len(replayed["redirected_main_figure_links"]) == 7
    assert len(replayed["preserved_source_pages"]) == 28
    with fitz.open(tmp_path / "replayed/document.pdf") as after, fitz.open(
            RESULTS / "publication_main_with_figures_20261004/document.pdf") as before:
        assert len(before) == len(after) == 32
        assert before.get_toc() == after.get_toc()
        for original, copied in zip(before, after):
            assert original.rect == copied.rect
            assert original.get_text("words") == copied.get_text("words")
            assert original.get_pixmap().samples == copied.get_pixmap().samples
            assert [(link["kind"], link.get("page"), link.get("file"), link.get("uri"), link["from"])
                    for link in original.get_links()] == [
                        (link["kind"], link.get("page"), link.get("file"), link.get("uri"), link["from"])
                        for link in copied.get_links()]
