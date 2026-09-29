import json
from pathlib import Path
import subprocess

import fitz
import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.print_manuscript_review import print_review


def setup(tmp_path):
    html = tmp_path / "source.html"
    html.write_text("<p>review</p>")
    browser = tmp_path / "browser"
    browser.write_text("fixture launcher")
    assets = tmp_path / "assets.json"
    assets.write_text(json.dumps(dict(status="manuscript_review_rendered",
                                     html=record(html), sources=[], targets=[])))
    return assets, tmp_path / "attempt", browser, html


def fake_browser(command, **kwargs):
    pdf = Path(next(arg.split("=", 1)[1] for arg in command if arg.startswith("--print-to-pdf=")))
    with fitz.open() as document:
        document.new_page().insert_text((50, 50), "review")
        document.save(pdf)
    return subprocess.CompletedProcess(command, 0, "printed", "")


def test_success_and_no_overwrite(tmp_path, monkeypatch):
    args = setup(tmp_path)
    monkeypatch.setattr(subprocess, "run", fake_browser)
    result = print_review(*args[:3])
    assert result["status"] == "verified_html_printed"
    assert result["page_count"] == 1
    assert result["visual_review_complete"] is False
    assert "--no-sandbox" not in result["command"]
    before = (args[1] / "document.pdf").read_bytes()
    with pytest.raises(FileExistsError):
        print_review(*args[:3])
    assert (args[1] / "document.pdf").read_bytes() == before


@pytest.mark.parametrize("kind", ["html", "receipt", "missing", "symlink"])
def test_reject_before_launch(tmp_path, monkeypatch, kind):
    assets, output, browser, html = setup(tmp_path)
    if kind == "html":
        html.write_text("changed")
    elif kind == "receipt":
        assets.write_text('{"status":"failed"}')
    elif kind == "missing":
        assets.unlink()
    else:
        output.symlink_to(tmp_path / "absent", target_is_directory=True)
    def forbidden(*args, **kwargs):
        pytest.fail("browser launched")
    monkeypatch.setattr(subprocess, "run", forbidden)
    with pytest.raises((ValueError, FileNotFoundError, FileExistsError)):
        print_review(assets, output, browser)
    assert not output.exists()


@pytest.mark.parametrize("kind", ["exit", "missing_pdf", "drift", "timeout"])
def test_failed_attempt_retained(tmp_path, monkeypatch, kind):
    assets, output, browser, html = setup(tmp_path)
    def failure(command, **kwargs):
        if kind == "timeout":
            raise subprocess.TimeoutExpired(command, 120)
        if kind == "drift":
            fake_browser(command, **kwargs)
            html.write_text("changed during print")
        return subprocess.CompletedProcess(command, 1 if kind == "exit" else 0, "", "diagnostic")
    monkeypatch.setattr(subprocess, "run", failure)
    with pytest.raises(Exception):
        print_review(assets, output, browser)
    result = json.loads((output / "print.json").read_text())
    assert result["status"] == "print_failed"
    assert "pdf" not in result
    assert result["publication_ready"] is False
