import json
from pathlib import Path

import pytest

from benchmark_tools import export_swiss_descriptive_strata as module

ROOT = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
COUNTS = ROOT / "qfo_fastoma_swiss_uncertainty_22098.json"
STRATA = ROOT / "corrected_swiss_sequence_strata_20260918.json"


@pytest.fixture
def inputs():
    return json.loads(COUNTS.read_text()), json.loads(STRATA.read_text())


def test_all_methods_bins_and_missingness(inputs):
    report, strata = inputs
    rows = module.build(report, strata)
    assert len(rows) == 8 * 11
    for row in rows:
        if row["method"] == "orthomcl_1_4" or row["families"] == 0:
            assert row["F1"] is None
        else:
            assert 0 <= row["F1"] <= 1
        if row["stratum"] == "all" and row["F1"] is not None:
            assert row["F1"] == pytest.approx(report["point_estimates"][row["method"]]["F1"], abs=1e-12)


def test_harmonic_macro_not_average_family_f1():
    values = [(.9, .1), (.1, .9)]
    assert module.scores(values)["F1"] == .5
    assert sum(module.scores([v])["F1"] for v in values)/2 == pytest.approx(.18)


def test_primary_bins_match_prior_admitted_analysis(inputs):
    previous = json.loads((ROOT / "qfo_corrected_swiss_primary_strata_21981.json").read_text())
    names = {"orthohmm_high_sensitivity": "high_sensitivity",
             "orthohmm_phylogeny_satellite_v2": "phylogenetic",
             "orthofinder_3_1_5_full": "orthofinder_full"}
    for row in module.build(*inputs):
        if row["method"] in names and row["stratum"] in ("lower_entropy", "higher_entropy"):
            prior = previous["bins"][row["stratum"].split("_")[0]]["point_estimates"][names[row["method"]]]
            for metric in module.METRICS:
                assert row[metric] == pytest.approx(prior[metric], abs=1e-12, rel=0)


@pytest.mark.parametrize("fault", ["count", "family", "bin", "aggregate", "nan", "missing", "duplicate"])
def test_invalid_input_rejected(inputs, fault):
    report, strata = inputs
    row = report["reconstructed_counts"]["methods"][0]
    if fault == "count": row["families"][0]["counts_without_prior"]["TP"] = -1
    elif fault == "family": row["families"].pop()
    elif fault == "bin": strata["primary_strata"]["lower_entropy"].append("unknown")
    elif fault == "aggregate": report["point_estimates"][row["method"]]["F1"] = 0.
    elif fault == "nan": row["families"][0]["statistics_with_prior"]["F1"] = float("nan")
    elif fault == "missing": report["point_estimates"]["orthomcl_1_4"] = {}
    elif fault == "duplicate": report["reconstructed_counts"]["methods"][-1] = row
    with pytest.raises(ValueError): module.build(report, strata)


def test_frozen_export_and_no_overwrite(tmp_path):
    output = tmp_path / "table"
    result = module.export(COUNTS, STRATA, output)
    assert result["new_inferential_claims"] is False
    assert len(result["outputs"]) == 2
    assert "NA" in (output / "scores.md").read_text()
    raw = (output / "scores.tsv").read_bytes()
    assert b"\r" not in raw and b"\t\n" not in raw
    with pytest.raises(FileExistsError): module.export(COUNTS, STRATA, output)
