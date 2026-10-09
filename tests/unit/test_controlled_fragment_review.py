"""Reuse the checked profile style without changing fragment evidence or HTML."""

import json
from pathlib import Path

import pytest

from benchmark_tools import render_controlled_fragment_review as current


def original(tmp_path):
    base = tmp_path / "benchmark_tools/results"
    base.mkdir(parents=True)
    heads = ("Suite", "Bin", "Families", "F1 Difference", "PPV Difference", "TPR Difference")
    table = "<table><tr>" + "".join("<th>" + s + "</th>" for s in heads) + "</tr>"
    table += "<tr>" + "<td>0</td>"*6 + "</tr>"
    table += ("<tr>" + "<td>0</td>"*6 + "</tr>")*22 + "</table>"
    text = '<html><head></head><body><h3 id="' + current.profile.ANCHOR + '">Profile</h3><p>Fixed bins</p>' + table
    text += '<h3 id="controlled-fragment-accuracy">Fragment</h3><p>Unchanged</p><table><tr><td>96.9773</td></tr></table></body></html>'
    html = base / "original.html"
    html.write_text(text)
    assets = base / "original.json"
    assets.write_text(json.dumps(dict(status="manuscript_review_rendered", publication_ready=False,
        html=current.record(html), sources=[], targets=[], render_command=["pandoc", "original.md"])))
    return base, html, assets, text


def test_scoped_recovery_inherits_actual_receipt_and_preserves_all_html(tmp_path):
    base, html, original_assets, text = original(tmp_path)
    output, assets = base / "corrected.html", base / "corrected.json"
    result = current.run(tmp_path, original_assets, output, assets)
    assert output.read_text().replace(current.profile.STYLE, "", 1) == text == html.read_text()
    assert result["render_command"] == ["pandoc", "original.md"]
    assert result["rendering_recovery"]["fragment_tables_changed"] is False
    assert result["rendering_recovery"]["scientific_html_unchanged"] is True
    assert result["rendering_recovery"]["all_columns_required"] == 6
    assert result["rendering_recovery"]["original_assets"] == current.record(original_assets)
    assert result["publication_ready"] is False


@pytest.mark.parametrize("change", ["modified_html", "changed_dom", "scope", "duplicate_style"])
def test_changed_or_invalid_review_refused(tmp_path, change):
    base, html, assets, text = original(tmp_path)
    receipt = json.loads(assets.read_text())
    if change == "modified_html": html.write_text(text + "changed")
    elif change == "scope": receipt["publication_ready"] = True
    else:
        html.write_text(text.replace("<td>0</td>", "", 1) if change == "changed_dom" else current.profile.styled(text))
        receipt["html"] = current.record(html)
    assets.write_text(json.dumps(receipt))
    with pytest.raises(ValueError): current.run(tmp_path, assets, base / "corrected.html", base / "corrected.json")
    assert not (base / "corrected.html").exists()


def test_existing_destination_guard_precedes_loading(tmp_path):
    base = tmp_path / "benchmark_tools/results"
    base.mkdir(parents=True)
    path = base / "exists.html"
    path.write_text("retained")
    with pytest.raises(FileExistsError): current.run(tmp_path, base / "absent.json", path, base / "new.json")
