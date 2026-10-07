"""Bounded four-cell manuscript generation and preservation of prior science."""

from copy import deepcopy
import json
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import prepare_four_cell_main_text as main
from benchmark_tools.check_native_qfo_supplement_presentation import table_rows

RESULTS = Path(main.__file__).parent / "results"


def inputs():
    return main.inputs(RESULTS)


def test_new_section_parses_and_all_tables_match_machine_readable_evidence():
    docs, refs = inputs()
    section = main.section(docs, refs)
    ast = json.loads(subprocess.check_output(["pandoc", "--from=markdown", "--to=json"], input=section, text=True))
    tables = [table_rows(block) for block in ast["blocks"] if block["t"] == "Table"]
    assert [len(rows) for rows in tables] == [6, 4, 9]
    cells = {r["cell"]: r for r in docs["snapshot"]["rows"]}
    for endpoint, row in zip(main.ENDPOINTS, tables[0]):
        assert row[0] == endpoint
        assert row[1] == ("F1" if endpoint in main.ENDPOINTS[:3] else "Sample mean" if endpoint == "FAS" else "Similarity")
        assert row[2:] == [f"{cells[c]['scores'][endpoint]:.8f}" for c in main.CELLS]
    for cell, row in zip(main.CELLS, tables[1]):
        assert row[0] == cell.upper().replace("_", "/")
        assert row[1] == f"{cells[cell]['submitted_pairs']:,}"
        assert row[2] == f"{100*cells[cell]['relation_coverage']:.4f}%"
    effects = {r["name"]: r for r in docs["binding"]["contrasts"]}
    for row, (name, metric) in zip(tables[2], [(c, m) for c in main.CONTRASTS for m in ("F1", "PPV", "TPR")]):
        value = effects[name]["metrics"][metric]
        assert row[2] == f"{100*value['difference']:+.4f}"
        assert row[3] == "[" + ", ".join(f"{100*x:.4f}" for x in value["bonferroni_percentile_ci"]) + "]"
        assert row[4] == "/".join(str(value[k]) for k in ("family_wins", "family_ties", "family_losses"))
    assert "three remain\nunavailable" in section and "other 11 remain unavailable" in section
    assert "not a new duplication exclusion" in section and "not biological correctness" in section


def test_all_other_scientific_body_is_unchanged_and_parent_bytes_retained():
    docs, refs = inputs()
    original = docs["parent"][docs["parent"].index("## Abstract\n"):]
    text = main.manuscript(docs, refs)
    restored = text[text.index("## Abstract\n"):]
    restored = restored.replace(main.section(docs, refs), original[original.index(main.START):original.index(main.END)], 1)
    restored = restored.replace("For P0/C0, every R-on GO scored pair", "Every R-on GO scored pair", 1)
    restored = restored.replace("They do not contain this later four-cell native QfO source snapshot.",
                                "They do not contain this later three-cell/complete-path\nsource snapshot.", 1)
    restored = restored.replace("../prepare_four_cell_main_text.py", "../prepare_error_strata_main_text.py", 1)
    assert restored == original
    assert main.record(RESULTS / main.PINS["parent"][0]) == refs["parent"]
    assert "Not submission-ready" in text and "prior P0 strata" in text


@pytest.mark.parametrize("change", ["schema", "ready", "missing_as_zero", "extra_admission", "duplicate_cell",
    "score", "mean", "denominator", "coverage", "pairs", "failed_timing", "shared_timing", "binding",
    "figure_reader", "draws", "multiplicity", "extra_contrast", "profile_readback", "localization",
    "tree_count", "tree_species"])
def test_inconsistent_or_inflated_evidence_refused(change):
    docs, refs = inputs()
    if change == "schema": docs["snapshot"]["schema"] = "previous"
    elif change == "ready": docs["rational"]["publication_ready"] = True
    elif change == "missing_as_zero": docs["snapshot"]["rows"][5]["scores"]["FAS"] = 0
    elif change == "extra_admission": docs["snapshot"]["rows"][5]["accuracy_admitted"] = True
    elif change == "duplicate_cell": docs["snapshot"]["rows"][5]["cell"] = main.CELLS[0]
    elif change == "score": docs["snapshot"]["rows"][4]["scores"]["FAS"] = float("nan")
    elif change == "mean": docs["snapshot"]["rows"][4]["secondary_mean"] = 0
    elif change == "denominator": docs["snapshot"]["rows"][4]["input_accessions"] = 10
    elif change == "coverage": docs["snapshot"]["rows"][4]["relation_coverage"] = .9
    elif change == "pairs": docs["snapshot"]["rows"][4]["submitted_pairs"] = True
    elif change == "failed_timing": docs["snapshot"]["rows"][1]["timing_eligible"] = True
    elif change == "shared_timing": docs["snapshot"]["rows"][4]["scientific_timings_admitted"] = True
    elif change == "binding": docs["binding"]["snapshot"] = {}
    elif change == "figure_reader": docs["figure_reader"]["exact_table_readback"] = False
    elif change == "draws": docs["binding"]["new_bootstrap_draws"] = True
    elif change == "multiplicity": docs["binding"]["multiplicity_endpoints"] = 9
    elif change == "extra_contrast": docs["binding"]["contrasts"][-1]["metrics"] = {}
    elif change == "profile_readback": docs["rational"]["contrast"] = {}
    elif change == "localization": docs["localization"]["summary"] = {"positive_paralogy_exclusion": 4}
    elif change == "tree_count": docs["tree_reader"]["tree_leaves_checked"] = 69
    else: docs["tree_reader"]["species_tree_bytes_identical"] = True
    with pytest.raises(ValueError):
        main.manuscript(docs, refs)


@pytest.mark.parametrize("change", ["duplicate_start", "missing_end", "wrong_ancillary"])
def test_parent_anchor_change_refused(change):
    docs, refs = inputs()
    if change == "duplicate_start": docs["parent"] += main.START
    elif change == "missing_end": docs["parent"] = docs["parent"].replace(main.END, "changed\n")
    else: docs["parent"] = docs["parent"].replace("Every R-on GO scored pair", "Every other scored pair")
    with pytest.raises(ValueError):
        main.manuscript(docs, refs)


def test_existing_outputs_refused_without_scanning_evidence(tmp_path):
    output, receipt = tmp_path / "main.md", tmp_path / "receipt.json"
    output.write_text("retain\n")
    with pytest.raises(ValueError, match="fresh"):
        main.generate(tmp_path, output, receipt)
    assert output.read_text() == "retain\n" and not receipt.exists()


def test_mismatched_relative_link_directory_refused(tmp_path):
    with pytest.raises(ValueError, match="beside"):
        main.generate(RESULTS, tmp_path / "main.md", tmp_path / "receipt.json")


def test_changed_input_digest_refused(tmp_path):
    (tmp_path / main.PINS["parent"][0]).write_text("not the preserved parent\n")
    with pytest.raises(ValueError, match="Changed manuscript evidence"):
        main.inputs(tmp_path)
