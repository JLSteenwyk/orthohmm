from pathlib import Path
from collections import Counter

import numpy as np
import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_canonical_downstream import filename_block_order, ordered_checkpoint, commands
from benchmark_tools.run_membership_hillclimb import argument
from orthohmm.accuracy import write_accuracy_checkpoint, load_accuracy_checkpoint


def test_filename_order_is_not_numeric_species_order():
    species = np.array([0, 0, 1, 1], dtype=np.int32)
    q = np.array([0, 1, 2, 3], dtype=np.int32)
    t = np.array([2, 3, 0, 1], dtype=np.int32)
    starts, ends, order = filename_block_order(species, q, t, {0: "z.fa", 1: "a.fa"})
    assert starts.tolist() == [0, 2] and ends.tolist() == [2, 4]
    assert order == [1, 0]


@pytest.mark.parametrize("classes", [{0: "a"}, {0: "a", 1: "a"}, {0: "a", 1: "b", 2: "c"}])
def test_filename_order_rejects_nonbijective_classes(classes):
    with pytest.raises(ValueError, match="bijective"):
        filename_block_order(np.array([0, 1]), np.array([0, 1]), np.array([1, 0]), classes)


def test_filename_order_rejects_interleaved_blocks():
    with pytest.raises(ValueError, match="contiguous"):
        filename_block_order(np.array([0, 1]), np.array([0, 1, 0]), np.array([1, 0, 1]), {0: "a", 1: "b"})


def test_checkpoint_order_preserves_duplicates_scores_and_internal_rows(tmp_path):
    names = ["a", "b", "c", "d"]
    species = np.array([0, 0, 1, 1], dtype=np.int32)
    q = np.array([0, 0, 1, 2, 3], dtype=np.int32)
    t = np.array([2, 2, 3, 0, 1], dtype=np.int32)
    scores = np.array([1.25, 1.25, 7.1, 3.5, 2.2])
    source = write_accuracy_checkpoint(str(tmp_path / "source"), names, species, q, t, scores)
    original = record(source / "manifest.json")
    owners = {"a": "z.fa", "b": "z.fa", "c": "a.fa", "d": "a.fa"}
    result = ordered_checkpoint(source, original["sha256"], owners, tmp_path / "output")
    loaded_names, loaded_species, nq, nt, ns = load_accuracy_checkpoint(result["checkpoint"])
    assert loaded_names == names and np.array_equal(loaded_species, species)
    assert nq.tolist() == [2, 3, 0, 0, 1]
    assert Counter(zip(q, t, scores)) == Counter(zip(nq, nt, ns))
    assert ns.tobytes() == scores[[3, 4, 0, 1, 2]].tobytes()
    assert record(source / "manifest.json") == original
    assert result["species_pair_blocks"] == 2


def test_checkpoint_requires_untampered_manifest_and_full_ownership(tmp_path):
    source = write_accuracy_checkpoint(str(tmp_path / "source"), ["a", "b"], [0, 1], [0, 1], [1, 0], [1., 1.])
    sha = record(source / "manifest.json")["sha256"]
    with pytest.raises(ValueError):
        ordered_checkpoint(source, "0" * 64, {"a": "z", "b": "a"}, tmp_path / "bad")
    with pytest.raises(ValueError, match="universes"):
        ordered_checkpoint(source, sha, {"a": "z"}, tmp_path / "missing")
    assert not (tmp_path / "bad").exists() and not (tmp_path / "missing").exists()


def test_commands_preserve_scientific_settings_and_use_own_checkpoint():
    graph, native = commands(Path("/frozen"), {"cpu": 16},
        {"fasta_directory": "/inputs", "aligner": "mafft", "tree_builder": "FastTree", "baseline_output": "/baseline"},
        Path("/new/canonical"), Path("/new/canonical/ordered/cp"), "abc")
    assert argument(graph, "--accuracy-checkpoint") == "/new/canonical/ordered/cp"
    assert argument(graph, "--checkpoint-sha256") == "abc"
    assert argument(graph, "--profile-iterations") == "1"
    assert argument(native, "--species-tree-mode") == "infer"
    assert argument(native, "--gene-tree-support-mode") == "legacy"
    assert argument(native, "--checkpoint-source") == "/baseline"
    assert argument(native, "--membership-constraints").startswith("/new/canonical/")
    assert "--species-tree" not in native and "--official-benchmark" not in graph


def test_failed_controls_skip_variants_but_continue_other_datasets(tmp_path, monkeypatch):
    import benchmark_tools.run_canonical_downstream as module
    from benchmark_tools.run_membership_hillclimb import write_json

    source = Path(module.__file__).resolve().parent.parent
    output = tmp_path / "output"
    output.mkdir()
    fasta = tmp_path / "a.fa"
    fasta.write_text(">gene\nACDEFG\n")
    cases = [{"dataset": dataset, "fastas": [record(fasta)], "checkpoint": "/unused",
              "checkpoint_manifest": {"sha256": "unused"}, "fasta_directory": str(tmp_path),
              "aligner": "mafft", "tree_builder": "FastTree", "baseline_output": "/unused"}
             for dataset in module.DATASETS]
    plan = output / "execution_plan.json"
    write_json(plan, {"repo": str(tmp_path), "source_root": str(source), "output": str(output),
                     "cpu": 1, "sources": [], "inputs": [], "runtime": {"files": []},
                     "tool_path": "/unused", "source_revision": "test", "scope": "mocked", "cases": cases})
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(module, "runtime", lambda: {"files": []})
    monkeypatch.setattr(module, "host_activity", lambda: "mocked host")
    monkeypatch.setattr(module, "GENES", dict.fromkeys(module.DATASETS, 1))
    monkeypatch.setattr(module, "audit", lambda *_: {})
    monkeypatch.setattr(module, "load_accuracy_checkpoint", lambda *a, **k: (
        ["gene"], np.array([0], dtype=np.int32), np.array([0], dtype=np.int32),
        np.array([0], dtype=np.int32), np.array([1.])))
    executed = []

    def failed_stage(command, folder, *args):
        executed.append(folder)
        return {"returncode": 1}

    def forbidden(*args, **kwargs):
        raise AssertionError("Failed control must not order variants or score")

    monkeypatch.setattr(module, "timed_stage", failed_stage)
    monkeypatch.setattr(module, "ordered_checkpoint", forbidden)
    monkeypatch.setattr(module, "score_orthobench", forbidden)
    result = module.run(plan)
    assert result["status"] == "finished_with_explicit_gaps"
    assert len(executed) == 3
    assert [row["status"] for row in result["arms"]] == ["failed", "not_run_failed_control_gate"] * 3
    assert all(row["scores"] == {} for row in result["arms"])
