import importlib
import inspect
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
NEW_COUNTS = ROOT / "qfo_recovered_swiss_uncertainty_22178.json"
NEW_SHA = "a599749d66433ec211ed1ab0e3a6a2a75eb99cb0dc57abc47c2d267930cbbe34"


@pytest.mark.parametrize("kind", ["descriptive", "identity", "fragment", "duplication"])
def test_explicit_counts_preserve_previous_rows_and_reject_wrong_hash(kind, tmp_path):
    module = importlib.import_module(f"benchmark_tools.export_swiss_{kind}_strata")
    old = json.loads((ROOT / f"swiss_{kind}_strata_20260923/manifest.json").read_text())
    feature = old["inputs"][1]
    args = [NEW_COUNTS, ROOT / Path(feature["path"]).name]
    if kind == "fragment":
        args.append(feature["sha256"])
    output = tmp_path / "new"
    assert inspect.signature(module.export).parameters["counts_sha"].default == module.COUNTS_SHA
    with pytest.raises(ValueError):
        module.export(*args, output, counts_sha="0" * 64)
    assert not output.exists()
    result = module.export(*args, output, counts_sha=NEW_SHA)
    previous = {(r["method"], r["stratum"]): r for r in old["rows"]}
    for row in result["rows"]:
        if row["method"] != "orthomcl_1_4":
            assert row == previous[row["method"], row["stratum"]]
        elif row["families"]:
            assert row["status"] == "descriptive"
            assert 0 <= row["F1"] <= 1
        else:
            assert row["F1"] is None
    assert result["new_inferential_claims"] is False
