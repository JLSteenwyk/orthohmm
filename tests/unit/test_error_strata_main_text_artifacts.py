"""Actual source snapshots and an explicitly reversible editorial-only revision."""

import hashlib
import json
from pathlib import Path

from benchmark_tools import prepare_error_strata_main_text as generator

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "benchmark_tools/results"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_actual_generated_source_and_bibliography_match_prospective_generator():
    receipt = json.loads((BASE / "error_strata_main_generation_20261007_v1.json").read_text())
    docs = {k: json.loads((BASE / name).read_text()) for k, (name, _) in generator.PINS.items()
            if name.endswith(".json")}
    parent = (BASE / generator.PINS["parent"][0]).read_text()
    expected, blocks = generator.extend(parent, docs)
    path = BASE / "PUBLICATION_MAIN_TEXT_20261007_v1.md"
    assert path.read_text() == expected
    assert digest(path) == receipt["output"]["sha256"]
    assert path.stat().st_size == receipt["output"]["bytes"]
    bibliography = json.loads((BASE / "publication_bibliography_20261007_v1.csl.json").read_text())
    assert bibliography[:-1] == docs["bibliography"]
    assert bibliography[-1] == generator.citation_entry(docs["citation"])
    assert receipt["blocks_sha256"] == {k: hashlib.sha256(v.encode()).hexdigest() for k, v in blocks.items()}
    assert receipt["bibliography_entries_added"] == 1
    assert receipt["source_commit"] == "522cc88239f4e1f08a7838d132058518b38b689d"
    assert digest(generator.__file__) == receipt["source"]["sha256"]


def test_actual_editorial_revision_is_only_declared_spacing_and_version_note():
    receipt = json.loads((BASE / "error_strata_main_editorial_20261007_v2.json").read_text())
    parent = ROOT / receipt["parent"]["path"]
    output = ROOT / receipt["output"]["path"]
    assert digest(parent) == receipt["parent"]["sha256"] and parent.stat().st_size == receipt["parent"]["bytes"]
    assert digest(output) == receipt["output"]["sha256"] and output.stat().st_size == receipt["output"]["bytes"]
    expected = parent.read_text()
    original = (BASE / generator.PINS["parent"][0]).read_text()
    assert len(receipt["edits"]) == 10
    for edit in [*receipt["edits"], receipt["version_note"]]:
        assert edit["occurrences"] == 1 and expected.count(edit["from"]) == 1
        if edit in receipt["edits"]:
            assert edit["from"] not in original
        expected = expected.replace(edit["from"], edit["to"], 1)
    assert output.read_text() == expected
    assert [l for l in parent.read_text().splitlines() if l.startswith("|")] == [
        l for l in output.read_text().splitlines() if l.startswith("|")]


def test_generation_and_editorial_receipts_do_not_claim_review_or_readiness():
    generation = json.loads((BASE / "error_strata_main_generation_20261007_v1.json").read_text())
    editorial = json.loads((BASE / "error_strata_main_editorial_20261007_v2.json").read_text())
    for receipt in (generation, editorial):
        assert receipt["publication_ready"] is False and receipt["new_render_or_package_proved"] is False
    assert generation["old_evidence_modified"] is False and generation["native_inference_reexecuted"] is False
    assert generation["new_bootstrap_draws"] == 0
    assert editorial["benchmark_or_scoring_reexecuted"] is False
    assert editorial["scientific_scope_unchanged"] is True and editorial["numerical_table_bytes_unchanged"] is True
    assert digest(BASE / generator.PINS["parent"][0]) == generator.PINS["parent"][1]


def test_new_citation_and_evidence_links_are_present_without_unreferenced_supplement_claims():
    draft = (BASE / "PUBLICATION_MAIN_TEXT_20261007_v2.md").read_text()
    for target in ("PUBLICATION_MAIN_TEXT_20261007_v1.md", "SWISS_MODEL_DIVERGENCE_PROTOCOL_20261007.md",
                   "swiss_model_divergence_strata_20261007_v1/TABLE.md", "swiss_model_divergence_evidence_23932_v1.tar.gz",
                   "native_qfo_three_cell_strata_figure_20261007_v1/fixed_stratum_contrasts.pdf"):
        assert target in draft and (BASE / target).is_file()
    assert "[@iqtree3_2026]" in draft and "IQ-TREE 3.0.1" in draft
    assert "all 11 nonempty" in draft and "66 points" in draft
    assert "all 10,765 distances" in draft
    assert "not calibrated biological time" in draft
