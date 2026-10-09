"""Verify scoped print recovery preserves every scientific HTML byte."""

import json
from pathlib import Path

import pytest

from benchmark_tools import render_profile_stratum_review as current


ROOT = Path(__file__).resolve().parents[2]


def invented():
    header = "".join("<th>" + s + "</th>" for s in
                     ("Suite", "Bin", "Families", "F1 Difference", "PPV Difference", "TPR Difference"))
    rows = "".join("<tr>" + "".join("<td>" + s + "</td>" for s in
                                 ("sequence", "invented_"+str(i), "0", "NA", "NA", "NA")) + "</tr>" for i in range(23))
    return '<html><head></head><body><h3 id="' + current.ANCHOR + '">Test</h3><p>Test-only</p><table><tr>' + header + "</tr>" + rows + "</table><table><tr><td>Unchanged other table</td></tr></table></body></html>"


def test_scoped_style_only_and_original_body_restorable():
    html = invented()
    updated = current.styled(html)
    assert updated.replace(current.STYLE, "", 1) == html
    assert "table-layout: fixed" in current.STYLE and "overflow-wrap: anywhere" in current.STYLE
    assert current.SELECTOR in current.STYLE and "nth-child(n+4)" in current.STYLE
    assert 'style id="profile-stratum-print-recovery"' in updated


@pytest.mark.parametrize("change", ["duplicate_anchor", "missing_column", "row_count", "wrong_selector", "already_styled", "head"])
def test_unsupported_dom_not_changed(change):
    html = invented()
    if change == "duplicate_anchor": html += '<h3 id="'+current.ANCHOR+'">Duplicate</h3>'
    elif change == "missing_column": html = html.replace("<th>TPR Difference</th>", "")
    elif change == "row_count": html = html.replace("</table>", "<tr><td>Extra</td></tr></table>", 1)
    elif change == "wrong_selector": html = html.replace("<p>Test-only</p>", "<div>Test-only</div>")
    elif change == "already_styled": html = current.styled(html)
    elif change == "head": html += "</head>"
    with pytest.raises(ValueError): current.styled(html)


def test_bound_recovery_and_changed_original_refusal(tmp_path, monkeypatch):
    base = tmp_path / "benchmark_tools/results"
    base.mkdir(parents=True)
    html = base / "invented.html"
    html.write_text(invented())
    original_assets = base / "original.json"
    original_assets.write_text(json.dumps(dict(status="manuscript_review_rendered", html=current.record(html),
        publication_ready=False, sources=[], targets=[], render_command=["invented"])))
    monkeypatch.setattr(current, "PINS", {"html":(html.name,current.record(html)["sha256"]),
                                          "assets":(original_assets.name,current.record(original_assets)["sha256"])})
    output, assets = base / "new.html", base / "new.json"
    result = current.run(tmp_path, output, assets)
    assert output.read_text().replace(current.STYLE, "", 1) == html.read_text()
    assert result["rendering_recovery"]["scientific_html_unchanged"] is True
    assert result["rendering_recovery"]["profile_rows"] == 23 and result["publication_ready"] is False
    html.write_text(html.read_text() + " ")
    with pytest.raises(ValueError): current.run(tmp_path, base / "refused.html", base / "refused.json")
    assert not (base / "refused.html").exists()


def test_overwrite_guard_before_loading(tmp_path):
    base = tmp_path / "benchmark_tools/results"
    occupied = base / "occupied"
    occupied.mkdir(parents=True)
    with pytest.raises(FileExistsError): current.run(tmp_path, occupied, base / "unused.json")


def test_actual_dom_structure_without_printing_or_source_mutation():
    path = ROOT / "benchmark_tools/results/native_qfo_four_cell_strata_manuscript_20261009_v1_review.html"
    original = path.read_text()
    assert current.styled(original).replace(current.STYLE, "", 1) == original
