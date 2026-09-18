from copy import deepcopy

import pytest

from benchmark_tools.prepare_qfo_factorial_pairs import reusable_stage, write_native_pairs


def reuse_fixture(profile=False):
    partition = {"path": "/candidate", "sha256": "partition", "bytes": 100}
    cell = {"reconciliation": False, "candidate_expansion": False, "profile_expansion": profile}
    labels = ["multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"]
    stages = [{"index": i, "stage": name, "partition": {**partition, "path": "/historical"},
               "retained_pairs": 8, "total_pairs": 10, "removed_mapping_pairs": 2} for i, name in enumerate(labels)]
    recovered = {"status": "four_stage_pairs_prepared_unscored", "accuracy_evaluated": False,
                 "input_fastas": [{"path": "/input"}], "mapping": {"sha256": "mapping"}, "stages": stages}
    return cell, {"candidate_partition": partition}, recovered, recovered["input_fastas"], recovered["mapping"]


@pytest.mark.parametrize("profile,index", [(False, 1), (True, 3)])
def test_reuses_only_identical_refined_partition(profile, index):
    assert reusable_stage(*reuse_fixture(profile))["index"] == index


@pytest.mark.parametrize("mutation", ["reconciliation", "expansion", "partition", "mapping", "fastas", "count", "stage", "scored"])
def test_invalid_reuse_rejected(mutation):
    cell, arm, recovered, fastas, mapping = deepcopy(reuse_fixture())
    if mutation == "reconciliation":
        cell["reconciliation"] = True
    elif mutation == "expansion":
        cell["candidate_expansion"] = True
    elif mutation == "partition":
        arm["candidate_partition"]["sha256"] = "changed"
    elif mutation == "mapping":
        mapping = {"sha256": "other"}
    elif mutation == "fastas":
        fastas = []
    elif mutation == "count":
        recovered["stages"][1]["removed_mapping_pairs"] = 9
    elif mutation == "stage":
        recovered["stages"][1]["stage"] = "multipass"
    else:
        recovered["accuracy_evaluated"] = True
    with pytest.raises(ValueError):
        reusable_stage(cell, arm, recovered, fastas, mapping)


def native(tmp_path, body):
    path = tmp_path / "native.tsv"
    path.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\n" + body)
    return path


def test_native_pairs_are_stripped_not_expanded(tmp_path):
    path = native(tmp_path, "sp|Z|X\tsp1\ttr|A|Y\tsp2\n")
    output = tmp_path / "pairs.tsv"
    assert write_native_pairs(path, output, {"sp|Z|X": "sp1", "tr|A|Y": "sp2", "tr|B|Q": "sp3"}, 1) == 1
    assert output.read_text() == "A\tZ\n"


@pytest.mark.parametrize("body,owners,count", [
    ("A\tx\tB\ty\nA\tx\tB\ty\n", {"A": "x", "B": "y"}, 2),
    ("B\ty\tA\tx\n", {"A": "x", "B": "y"}, 1),
    ("A\tx\tB\ty\n", {"A": "x", "B": "y"}, 2),
    ("A\tx\tB\ty\n", {"A": "x", "sp|A|X": "x", "B": "y"}, 1),
    ("A\tx\tB\tx\n", {"A": "x", "B": "x"}, 1)])
def test_native_conversion_rejects_invalid_input(tmp_path, body, owners, count):
    with pytest.raises(ValueError):
        write_native_pairs(native(tmp_path, body), tmp_path / "pairs.tsv", owners, count)
