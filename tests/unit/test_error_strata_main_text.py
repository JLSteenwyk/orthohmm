"""Bounded manuscript additions: preserve every pre-existing scientific byte."""

import copy
import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools import prepare_error_strata_main_text as generator

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "benchmark_tools/results"


@pytest.fixture
def inputs():
    docs = {k: json.loads((BASE / name).read_text()) for k, (name, _) in generator.PINS.items()
            if name.endswith(".json")}
    parent = (BASE / generator.PINS["parent"][0]).read_text()
    return parent, docs


def test_all_current_direct_input_pins_match():
    for name, sha in generator.PINS.values():
        assert hashlib.sha256((BASE / name).read_bytes()).hexdigest() == sha


def test_removing_only_new_blocks_and_selector_restores_parent_exactly(inputs):
    parent, docs = inputs
    revised, blocks = generator.extend(parent, docs)
    assert revised != parent and revised.count(generator.HEADER) == 1
    restored = revised.replace(generator.HEADER, "", 1)
    for text in blocks.values():
        assert restored.count(text) == 1
        restored = restored.replace(text, "", 1)
    restored = restored.replace(generator.NEW_SELECTOR, generator.OLD_SELECTOR, 1)
    assert restored == parent
    assert revised.startswith(parent.split("\n\n", 1)[0])


def test_all_six_distance_contrasts_all_metrics_and_sizes_in_manuscript(inputs):
    parent, docs = inputs
    text, _ = generator.extend(parent, docs)
    for row in docs["distance"]["differences"]:
        values = " | ".join(f"{row[m + '_pp']:+.3f}" for m in ("F1", "PPV", "TPR"))
        assert text.count("| " + values + " |") == 1
    assert text.count("| Lower/equal distance | 9 |") == 2
    assert text.count("| Higher distance | 9 |") == 2
    for statement in ("563canonical", "all10765distances", "all11nonempty", "66points",
                      "five empty bins", "four redundant all-family copies", "failed"):
        assert statement in text


@pytest.mark.parametrize("change", ["family", "count", "partial", "model", "timing", "scope", "figure", "sign"])
def test_inconsistent_or_inflated_inherited_science_rejected(inputs, change):
    _, docs = inputs
    docs = copy.deepcopy(docs)
    if change == "family":
        docs["distance"]["memberships"]["APP"] = []
    elif change == "count":
        docs["distance"]["family_rows"][0]["counts_without_prior"]["TP"] += 1
    elif change == "partial":
        docs["distance_reader"]["differences_checked"] = 5
    elif change == "model":
        docs["distance"]["model"] = "LG+G4"
    elif change == "timing":
        docs["distance"]["cells"][1]["timing_eligible"] = True
    elif change == "scope":
        docs["feature_reader"]["independent_confirmation"] = True
    elif change == "figure":
        docs["figure_review"]["full_pdf_page_raster_actually_viewed"] = False
    elif change == "sign":
        row = next(r for r in docs["fixed"]["differences"] if r["contrast"] == "C_at_P0_R0"
                   and r["status"] == "descriptive")
        row["PPV"] = 1
    with pytest.raises(ValueError):
        generator.blocks(docs)


def test_uncertainty_calibration_exposure_and_runtime_scope_are_explicit(inputs):
    parent, docs = inputs
    text, _ = generator.extend(parent, docs)
    prose = " ".join(text.split())
    for phrase in ("not confidence intervals", "not calibrated biological time",
                   "known ancestral history", "not independent confirmation",
                   "no new intervals or significance claim", "not complete ancestral histories"):
        assert phrase in prose
    assert "does not establish a new HTML/PDF" in generator.HEADER
    assert "TotalCPU and MaxRSS were unavailable" in text
    assert "no software update or rerun follows" in text
    assert "not submission-ready" in text.lower()


def test_changed_anchor_and_double_integration_fail(inputs):
    parent, docs = inputs
    with pytest.raises(ValueError, match="anchor"):
        generator.extend(parent.replace(generator.ANCHORS["results"], "### Changed\n"), docs)
    revised, _ = generator.extend(parent, docs)
    with pytest.raises(ValueError):
        generator.extend(revised, docs)


def test_current_citation_is_exact_publisher_metadata_not_old_submitted_notice(inputs):
    _, docs = inputs
    entry = generator.citation_entry(docs["citation"])
    deposited = docs["citation"]["message"]
    assert entry["id"] == "iqtree3_2026"
    assert entry["DOI"] == "10.1093/molbev/msag117" and entry["issued"]["date-parts"][0][0] == 2026
    assert entry["author"] == [{k: a[k] for k in ("given", "family")} for a in deposited["author"]]
    assert len(entry["author"]) == 16
    assert entry["id"] not in {e["id"] for e in docs["bibliography"]}
    assert "Submitted" not in json.dumps(entry)
    entries = [*docs["bibliography"], entry]
    assert entries[:-1] == docs["bibliography"]


def test_existing_output_refused_before_inputs_or_git_read(tmp_path):
    with pytest.raises(ValueError, match="Existing generation"):
        generator.generate(tmp_path, tmp_path, tmp_path / "newbib", tmp_path / "receipt", "a" * 40)
