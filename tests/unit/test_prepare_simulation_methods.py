from pathlib import Path

import pytest

from benchmark_tools.prepare_simulation_methods import commands, prepare, reuse_comparators


def test_frozen_commands_and_diagnostic_parent():
    dataset = {"input": "/data/input", "truth": "/secret_reference/truth.json"}
    result = commands(dataset, Path("/out"), Path("/frozen"), Path("/python"), Path("/of"))
    assert set(result) == {"orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full", "orthofinder_sequence_only"}
    for method in ("orthohmm_high_sensitivity", "orthohmm_satellite_v2"):
        argv = result[method]["argv"]
        for flag, value in (("--cpu", "4"), ("--threads-per-worker", "4"),
                            ("--accuracy-profile", "high_sensitivity"), ("--cpm-resolution", "0.1")):
            assert argv[argv.index(flag) + 1] == value
        assert "/secret_reference/truth.json" not in argv
    high = result["orthohmm_high_sensitivity"]["argv"]
    phylo = result["orthohmm_satellite_v2"]["argv"]
    assert "--phylogeny" not in high
    assert phylo[phylo.index("--species-tree-mode") + 1] == "infer"
    assert phylo[phylo.index("--species-tree-rooting") + 1] == "min_variance"
    assert phylo[phylo.index("--phylogeny-candidates") + 1] == "satellite_v2"
    assert result["orthofinder_full"]["argv"] == ["/of", "-f", "/out/orthofinder_full/input", "-t", "4", "-a", "4", "-S", "diamond"]
    assert result["orthofinder_full"]["copy_inputs_from"] == "/data/input"
    assert result["orthofinder_sequence_only"]["parent_method"] == "orthofinder_full"
    assert "argv" not in result["orthofinder_sequence_only"]


def test_wrong_generation_manifest_fails_before_external_inspection(tmp_path):
    path = tmp_path / "generation.json"
    path.write_text("{}")
    with pytest.raises(ValueError, match="Wrong frozen"):
        prepare(path, tmp_path, tmp_path / "python", tmp_path / "of", tmp_path / "out", tmp_path / "manifest.json")
    assert not (tmp_path / "manifest.json").exists()


def test_explicit_generation_hash_still_requires_full_panel(tmp_path):
    import hashlib
    path = tmp_path / "generation.json"
    path.write_text("{}")
    with pytest.raises(ValueError, match="complete frozen"):
        prepare(path, tmp_path, tmp_path / "python", tmp_path / "of", tmp_path / "out", tmp_path / "manifest.json",
                hashlib.sha256(path.read_bytes()).hexdigest())


def reuse_fixture():
    data = {"input": "/data/input", "truth": "/data/truth", "label": "baseline", "seed": 1}
    def manifest(output, root):
        return {"generation_manifest": {"sha256": "same"}, "core_commit": "frozen",
                "core_root": root, "tool_entrypoints": {"orthohmm_python": {"absolute_path": "/python"},
                                                          "orthofinder": {"absolute_path": "/of"}},
                "orthofinder_distribution": [], "environments": {},
                "datasets": [{**data, "methods": commands(data, Path(output), Path(root), Path("/python"), Path("/of"))}]}
    return manifest("/new", "/native"), manifest("/old", "/source")


def test_reuse_retains_old_comparator_paths_and_new_hmm_paths():
    report, previous = reuse_fixture()
    reuse_comparators(report, previous)
    assert report["execution_order"] == ["orthohmm_high_sensitivity", "orthohmm_satellite_v2"]
    methods = report["datasets"][0]["methods"]
    assert methods["orthofinder_full"]["output"] == "/old/orthofinder_full"
    assert methods["orthohmm_high_sensitivity"]["output"] == "/new/orthohmm_high_sensitivity"


@pytest.mark.parametrize("change", ["generation", "seed", "sensitivity", "comparator", "duplicate"])
def test_reuse_rejects_scientific_drift(change):
    report, previous = reuse_fixture()
    if change == "generation":
        previous["generation_manifest"]["sha256"] = "different"
    elif change == "seed":
        previous["datasets"][0]["seed"] = 2
    elif change == "sensitivity":
        previous["datasets"][0]["methods"]["orthohmm_high_sensitivity"]["argv"][-1] = "fast"
    elif change == "comparator":
        previous["datasets"][0]["methods"]["orthofinder_full"]["argv"][-1] = "blast"
    else:
        previous["datasets"].append(previous["datasets"][0])
    with pytest.raises(ValueError):
        reuse_comparators(report, previous)
