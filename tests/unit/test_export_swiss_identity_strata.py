import json
from pathlib import Path

import pytest

from benchmark_tools.export_swiss_identity_strata import rows_with_differences, REFERENCE

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def inputs():
    return (json.loads((RESULTS / "qfo_fastoma_swiss_uncertainty_22098.json").read_text()),
            json.loads((RESULTS / "corrected_swiss_identity_admission_22102.json").read_text()))


def test_retains_all_methods_bins_and_reference_differences():
    rows = rows_with_differences(*inputs())
    assert len(rows) == 32
    for row in rows:
        if row["stratum"] == "missing_identity" or row["method"] == "orthomcl_1_4":
            assert row["F1"] is row["delta_F1"] is None
        elif row["method"] == REFERENCE:
            assert row["delta_F1"] == row["delta_PPV"] == row["delta_TPR"] == 0
        else:
            ref = next(r for r in rows if r["method"] == REFERENCE and r["stratum"] == row["stratum"])
            assert row["delta_F1"] == row["F1"] - ref["F1"]


@pytest.mark.parametrize("problem", ["status", "scored", "duplicate", "missing", "extra_bin"])
def test_invalid_admission_rejected(problem):
    counts, admission = inputs()
    if problem == "status":
        admission["status"] = "prepared_only"
    elif problem == "scored":
        admission["prediction_statistics_evaluated"] = True
    elif problem == "duplicate":
        admission["strata"]["lower_identity"].append(admission["strata"]["higher_identity"][0])
    elif problem == "missing":
        admission["strata"]["lower_identity"].pop()
    else:
        admission["strata"]["extra"] = []
    with pytest.raises(ValueError):
        rows_with_differences(counts, admission)
