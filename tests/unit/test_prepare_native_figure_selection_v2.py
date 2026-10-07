import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import prepare_native_figure_selection_v2 as current


ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def inputs():
    return tuple(json.loads((ROOT / ref["path"]).read_text()) for ref in (current.PRIOR, current.PRINT))


def test_old_sixteen_figures_and_captions_are_preserved(inputs):
    prior, printed = inputs
    old = copy.deepcopy(prior)
    result = current.select(prior, printed)
    assert prior == old
    assert result["figures"][:16] == prior["figures"]
    assert len(result["figures"]) == 17 and result["main_pages"] == 19
    assert result["figures"][-1]["pdf"] == current.FIGURE
    assert result["publication_ready"] is False
    assert result["scientific_settings_or_results_changed"] is False
    assert "Four fresh cells and twelve contrasts remain unavailable" in result["figures"][-1]["caption"]


@pytest.mark.parametrize("change", ["pages", "pdf", "status", "ready", "order"])
def test_inconsistent_print_or_historical_inventory_refused(inputs, change):
    prior, printed = copy.deepcopy(inputs)
    if change == "pages": printed["page_count"] += 1
    elif change == "pdf": printed["pdf"]["sha256"] = "0" * 64
    elif change == "status": printed["status"] = "print_started"
    elif change == "ready": prior["publication_ready"] = True
    elif change == "order": prior["figures"].reverse()
    with pytest.raises(ValueError): current.select(prior, printed)


def test_selection_build_checks_all_actual_pdfs_and_refuses_overwrite(tmp_path, inputs):
    output = tmp_path / "selection.json"
    result = current.build(output)
    selection = json.loads(output.read_text())
    assert selection["figures"][:16] == inputs[0]["figures"]
    assert result["bytes"] == output.stat().st_size
    assert len(selection["provenance"]) == 5
    with pytest.raises(FileExistsError): current.build(output)


def test_different_input_bytes_cannot_pass_checksum(tmp_path):
    path = tmp_path / "input.json"
    path.write_text("{}")
    with pytest.raises(ValueError):
        current.checked(tmp_path, dict(path="input.json", bytes=2, sha256="0" * 64))
