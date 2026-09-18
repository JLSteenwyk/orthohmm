import json

import numpy as np
import pytest

from orthohmm.accuracy import write_accuracy_checkpoint
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.validate_high_sensitivity_outputs import validate, validate_metrics, validate_partition


def metrics():
    return {"status": "complete", "metadata": {"accuracy_profile": "high_sensitivity", "search_mode": "builtin",
            "clustering": "leiden", "cpm_resolution": 0.1, "substitution_matrix": "BLOSUM62",
            "evalue_threshold": 1e-4, "leiden_seed": 4, "cpu_budget": 32},
            "counts": {"genes": 2, "species": 2, "orthogroups": 1, "high_sensitivity_profiles": 1}}


@pytest.mark.parametrize("key,value", [("status", "running"), ("accuracy_profile", "standard"),
    ("cpu_budget", 8), ("leiden_seed", 9), ("evalue_threshold", 0.1)])
def test_reject_metrics(key, value):
    data = metrics()
    (data if key == "status" else data["metadata"])[key] = value
    with pytest.raises(ValueError):
        validate_metrics(data, 2, 2, 1)


@pytest.mark.parametrize("key", ["genes", "species", "orthogroups", "high_sensitivity_profiles"])
def test_reject_counts(key):
    data = metrics()
    data["counts"][key] = 0
    with pytest.raises(ValueError):
        validate_metrics(data, 2, 2, 1)


@pytest.mark.parametrize("group", ["OG0: a b\n", "OG0: a\n", "OG0: a a b\n", "OG0: a foreign\n"])
def test_real_checkpoint_and_groups(tmp_path, group):
    inputs = []
    for species, gene in (("s1", "a"), ("s2", "b")):
        path = tmp_path / f"{species}.fasta"
        path.write_text(f">{gene}\nACDE\n")
        inputs.append(path)
    output = tmp_path / "output"
    checkpoint = write_accuracy_checkpoint(str(output), ["a", "b"], np.array([0, 1], dtype=np.int32),
        np.array([0, 1], dtype=np.int32), np.array([1, 0], dtype=np.int32), np.array([10., 10.]))
    (output / "orthohmm_orthogroups.txt").write_text(group)
    (output / "orthohmm_working_res/orthohmm_edges_clustered.txt").write_text("a b\n")
    path = tmp_path / "metrics.json"
    data = metrics()
    data["metadata"].update(output_directory=str(output), fasta_directory=str(tmp_path))
    path.write_text(json.dumps(data))
    sha = record(checkpoint / "manifest.json")["sha256"]
    if group == "OG0: a b\n":
        result = validate(output, path, inputs, sha)
        assert result["genes"] == 2 and result["groups"] == 1
        assert result["numeric_checkpoint"]["summary"]["hits"] == 2
    else:
        with pytest.raises(ValueError):
            validate(output, path, inputs, sha)


@pytest.mark.parametrize("value", [None, True, 1.0, float("nan"), float("inf"), "1", -1])
def test_profile_count_requires_positive_integer(value):
    data = metrics()
    data["counts"]["high_sensitivity_profiles"] = value
    with pytest.raises(ValueError):
        validate_metrics(data, 2, 2, 1)


def test_actual_native_metrics_writer(tmp_path):
    from orthohmm.metrics import PipelineMetrics
    path = tmp_path / "native.json"
    fixture = metrics()
    with PipelineMetrics(str(path)) as writer:
        writer.add_metadata(**fixture["metadata"])
        writer.add_counts(**fixture["counts"])
    validate_metrics(json.loads(path.read_text()), 2, 2, 1)


def test_reject_flattened_metadata():
    data = metrics()
    data.update(data.pop("metadata"))
    with pytest.raises(ValueError, match="metadata"):
        validate_metrics(data, 2, 2, 1)


@pytest.mark.parametrize("raw,groups,valid", [
    ("a b\n", {"OG0": ["b", "a"], "OG1": ["c"]}, True),
    ("a b\nc\n", {"OG0": ["a", "b"], "OG1": ["c"]}, True),
    ("", {"OG0": ["a"], "OG1": ["b"], "OG2": ["c"]}, True),
    ("a b\n", {"OG0": ["a", "c"], "OG1": ["b"]}, False),
    ("a b\n", {"OG0": ["a", "b", "c"]}, False),
    ("a b\n", {"OG0": ["a"], "OG1": ["b"], "OG2": ["c"]}, False),
    ("a b\nb c\n", {"OG0": ["a", "b", "c"]}, False),
    ("a foreign\n", {"OG0": ["a", "b", "c"]}, False),
])
def test_cluster_export_equivalence(tmp_path, raw, groups, valid):
    path = tmp_path / "raw.txt"
    path.write_text(raw)
    if valid:
        result = validate_partition(path, groups, ["a", "b", "c"])
        assert result["raw_clusters"] + result["added_singletons"] == len(groups)
    else:
        with pytest.raises(ValueError):
            validate_partition(path, groups, ["a", "b", "c"])
