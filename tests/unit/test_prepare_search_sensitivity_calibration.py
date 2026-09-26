import json

import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.prepare_search_sensitivity_calibration import dataset, family_denominators, EVALUE_GRID


def test_directed_homology_includes_paralogs_but_not_same_species():
    assert family_denominators({"F": ["a", "b", "c"]}, {"a": "A", "b": "A", "c": "B"}) == {"F": 4}


def test_ineligible_family_retained_with_zero_denominator():
    assert family_denominators({"F": ["a"], "G": []}, {"a": "A"}) == {"F": 0, "G": 0}


@pytest.mark.parametrize("families", [{"F": ["a", "a", "b"]}, {"F": ["a"], "G": ["a", "b"]},
                                    {"F": ["a", "b", "x"]}, {"F": ["a"]}])
def test_invalid_family_membership(families):
    with pytest.raises(ValueError):
        family_denominators(families, {"a": "A", "b": "B"})


def test_fixed_grid():
    assert EVALUE_GRID == sorted(set(EVALUE_GRID))
    assert len(EVALUE_GRID) == 15
    assert EVALUE_GRID[0] == 1e-100 and EVALUE_GRID[-1] == 1


@pytest.fixture
def input_dataset(tmp_path):
    inputs = tmp_path / "input"
    inputs.mkdir()
    records = []
    for taxon, gene in (("A", "a"), ("B", "b")):
        path = inputs / f"{taxon}.fasta"
        path.write_text(f">{gene}\nACDE\n")
        item = record(path)
        item["path"] = f"input/{path.name}"
        records.append(item)
    truth = tmp_path / "truth.json"
    truth.write_text(json.dumps(dict(inputs=records, families={"F": ["a", "b"]},
                                    extant_genes=2, species=["A", "B"])))
    row = dict(condition="baseline", seed=20261101, input=str(inputs), truth=str(truth))
    return row, record(truth)["sha256"]


def test_dataset_inventory_and_split(input_dataset):
    row, sha = input_dataset
    result = dataset(row, sha)
    assert result["directed_homology_pairs"] == 2
    assert result["directed_nonhomology_pairs"] == 0
    assert result["split"] == "calibration"
    row["seed"] = 20261106
    assert dataset(row, sha)["split"] == "reporting"


@pytest.mark.parametrize("mutation", ["truth", "content", "extra", "foreign", "duplicate_input"])
def test_changed_dataset_rejected(input_dataset, mutation):
    from pathlib import Path
    row, sha = input_dataset
    truth = Path(row["truth"])
    if mutation == "truth":
        sha = "0" * 64
    elif mutation == "content":
        (Path(row["input"]) / "A.fasta").write_text(">a\nAAAA\n")
    elif mutation == "extra":
        (Path(row["input"]) / "extra.txt").write_text("extra")
    else:
        data = json.loads(truth.read_text())
        if mutation == "foreign":
            data["families"]["F"].append("x")
        else:
            data["inputs"].append(data["inputs"][0])
        truth.write_text(json.dumps(data))
        sha = record(truth)["sha256"]
    with pytest.raises(ValueError):
        dataset(row, sha)
