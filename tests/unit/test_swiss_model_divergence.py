"""New divergence arithmetic, safe attempt boundaries and native fixture smoke."""

import inspect
import math
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import prepare_swiss_model_divergence as writer
from benchmark_tools import readback_swiss_model_divergence as reader


@pytest.fixture
def tree(tmp_path):
    path = tmp_path / "tree.nwk"
    path.write_text("((A:0.1,B:0.2):0.3,C:0.4,D:0):999;\n")
    return path


def test_paths_and_independent_splits_agree_with_hand_arithmetic(tree):
    feature, pairs = writer.tree_features(tree, ["A", "B", "C", "D"])
    other, split_pairs = reader.edge_distances(tree, ["A", "B", "C", "D"])
    values = {("A", "B"): .3, ("A", "C"): .8, ("A", "D"): .4,
              ("B", "C"): .9, ("B", "D"): .5, ("C", "D"): .4}
    assert split_pairs == pytest.approx(values)
    assert {(p["gene_a"], p["gene_b"]): p["distance"] for p in pairs} == pytest.approx(values)
    assert feature["median_pair_distance"] == pytest.approx(.45)
    assert feature["tree_length"] == pytest.approx(1)
    for key in ("median_pair_distance", "mean_pair_distance", "tree_length"):
        assert feature[key] == pytest.approx(other[key])


@pytest.mark.parametrize("text", ["(A:-.1,B:.2);", "(A,B:.2);", "(A:.1,A:.2);",
                                      "(A:.1,X:.2);", "(A:.1,B:.2);(A:.3,B:.4);"])
@pytest.mark.parametrize("function", [writer.tree_features, reader.edge_distances])
def test_invalid_edges_tips_and_multiple_trees_rejected(tmp_path, text, function):
    path = tmp_path / "bad.nwk"
    path.write_text(text)
    with pytest.raises(ValueError):
        function(path, ["A", "B"])


def test_zero_edges_duplicate_sequences_and_missing_root_length_are_allowed(tmp_path):
    path = tmp_path / "zero.nwk"
    path.write_text("(A:0,B:0,C:.3);\n")
    feature, pairs = writer.tree_features(path, ["A", "B", "C"])
    assert feature["pairs"] == 3
    assert pairs[0]["distance"] == 0
    assert feature["median_pair_distance"] == .3
    assert reader.edge_distances(path, ["A", "B", "C"])[0] == feature


def test_even_median_and_ties_are_not_forcibly_balanced():
    features = {"A": {"median_pair_distance": 1}, "B": {"median_pair_distance": 1},
                "C": {"median_pair_distance": 1}, "D": {"median_pair_distance": 2}}
    cutoff, bins = writer.bins(features)
    assert cutoff == 1
    assert bins["lower_or_equal_median"] == ["A", "B", "C"]
    assert bins["higher_than_median"] == ["D"]
    cutoff, bins = writer.bins({k: {"median_pair_distance": 0} for k in features})
    assert cutoff == 0 and bins["higher_than_median"] == []


@pytest.mark.parametrize("value", [math.nan, math.inf, -1])
def test_invalid_bin_feature_rejected(value):
    with pytest.raises(ValueError):
        writer.bins({"A": {"median_pair_distance": value}})


@pytest.mark.parametrize("text", [">A\nAA-\n>A\nAA-\n", ">A\nAAA\n>B\nAA\n",
                                  ">A\nAA*\n>B\nAAA\n", ">A\n---\n>B\nAAA\n"])
def test_alignment_scope_alphabet_and_dimensions_guarded(tmp_path, text):
    path = tmp_path / "bad.faa"
    path.write_text(text)
    with pytest.raises(ValueError):
        writer.alignment_members(path, ["A", "B"], 3)


def test_pair_tsv_independent_numeric_missing_duplicate_and_tamper_checks(tmp_path):
    path = tmp_path / "pairs.tsv"
    expected = {("fam", "A", "B"): .3}
    path.write_text("family\tgene_a\tgene_b\tdistance\nfam\tA\tB\t0.3\n")
    assert reader.compare_pairs(path, expected) == 1
    for suffix in ("fam\tA\tB\t0.3\nfam\tA\tB\t0.3\n", "fam\tA\tB\t0.4\n",
                   "fam\tA\tB\tnan\n", ""):
        path.write_text("family\tgene_a\tgene_b\tdistance\n" + suffix)
        with pytest.raises(ValueError):
            reader.compare_pairs(path, expected)


def test_feature_stage_never_reads_prediction_scores():
    source = inspect.getsource(writer)
    assert "native_qfo_three_cell" not in source
    assert "counts_without_prior" not in source
    assert "Phylo.read" in source
    assert "tree.distance" in source
    independent = inspect.getsource(reader)
    assert "prepare_swiss_model_divergence" not in independent.split("def verify")[0]
    assert ".distance(" not in independent


def test_fixed_command_no_model_selection_support_resume_or_unbounded_threads():
    argv = writer.command("iqtree3", "retained.faa", "fresh/inference")
    assert argv == ["iqtree3", "-s", "retained.faa", "--seqtype", "AA", "-m", "WAG+G4",
                    "--seed", "20261007", "-T", "1", "--mem", "4G", "-keep-ident",
                    "--prefix", "fresh/inference"]
    assert not any(v in argv for v in ("MFP", "--redo", "-B", "--alrt", "AUTO"))


def test_existing_attempt_refuses_before_loading_inputs(tmp_path, monkeypatch):
    monkeypatch.setattr(writer, "load_inputs", lambda repo: pytest.fail("Must refuse existing attempt first"))
    with pytest.raises(ValueError, match="Existing attempt"):
        writer.execute(tmp_path, tmp_path, "a" * 40, 1)


def test_save_refuses_overwrite(tmp_path):
    path = tmp_path / "receipt.json"
    writer.save(path, {"status": "retained"})
    with pytest.raises(FileExistsError):
        writer.save(path, {"status": "replacement"})


def test_timeout_retained_and_subprocess_terminated(tmp_path):
    code, timeout = writer.run_family([sys.executable, "-c", "import time;time.sleep(60)"],
                                     tmp_path / "out", tmp_path / "err", timeout=.05)
    assert timeout is True and code != 0


def test_batch_syntax_and_scoped_resources():
    path = Path(writer.__file__).with_name("swiss_model_divergence_batch_20261007.sh")
    subprocess.run(["bash", "-n", str(path)], check=True)
    text = path.read_text()
    for scope in ("--cpus-per-task=2", "--mem=8G", "--time=04:00:00", "--no-requeue",
                  "OMP_NUM_THREADS=1", "-u LD_PRELOAD"):
        assert scope in text
    assert "exclusive" not in text and "quiet" not in text


@pytest.mark.skipif(not writer.BINARY.exists(), reason="Optional trusted local binary")
def test_one_failed_family_retained_no_subset_cutoff_and_other_families_continue(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    output = repo / "benchmarks/results/fresh"
    output.parent.mkdir(parents=True)
    members = {"A": ["gene1", "gene2"], "B": ["gene3", "gene4"]}
    runs = [dict(refog=f, alignment=dict(path="invented.faa"), columns=10) for f in members]
    monkeypatch.setattr(writer, "load_inputs", lambda repo: ({}, runs, members))
    monkeypatch.setattr(writer, "committed_sources", lambda repo, commit: [])
    for key, value in dict(SLURM_JOB_ID="1", SLURM_CPUS_PER_TASK="2", SLURM_MEM_PER_NODE="8192",
                           OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1").items():
        monkeypatch.setenv(key, value)
    calls = []

    def mocked_family(argv, stdout, stderr):
        family = Path(stdout).parent.name
        calls.append(family)
        Path(stdout).write_text("fixture\n")
        Path(stderr).write_text("")
        if family == "A":
            return 1, False
        Path(stdout).with_name("inference.iqtree").write_text("Model of substitution: WAG+G4\n")
        Path(stdout).with_name("inference.treefile").write_text("(gene3:.1,gene4:.2);\n")
        return 0, False

    monkeypatch.setattr(writer, "run_family", mocked_family)
    result = writer.execute(repo, output, "a" * 40, 1)
    assert calls == ["A", "B"]
    assert result["status"] == "selected_family_failure"
    assert result["failed_families"] == ["A"]
    assert result["median_family_distance"] is None and result["strata"] is None
    assert result["runs"][0]["exit_code"] == 1 and result["runs"][0]["features"] is None
    assert result["runs"][1]["status"] == "feature_constructed"
    assert (output / "A/result.json").exists() and (output / "report.json").exists()


@pytest.mark.skipif(not writer.BINARY.exists(), reason="Optional trusted local IQ-TREE smoke fixture")
def test_installed_iqtree_fixed_model_and_identical_tip_fixture(tmp_path):
    path = tmp_path / "invented.faa"
    base = "ACDEFGHIKLMNPQRSTVWY" * 4
    sequences = dict(A=base, B=base, C="C" + base[1:], D=base[:10] + "A" + base[11:],
                     E="ACCC" + base[4:])
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in sequences.items()))
    argv = writer.command(writer.BINARY, path, tmp_path / "inference")
    code, timeout = writer.run_family(argv, tmp_path / "out", tmp_path / "err", timeout=60)
    assert code == 0 and timeout is False
    report = (tmp_path / "inference.iqtree").read_text()
    assert "Model of substitution: WAG+G4" in report
    descriptors, _ = writer.tree_features(tmp_path / "inference.treefile", sorted(sequences))
    assert descriptors["pairs"] == 10
    other, _ = reader.edge_distances(tmp_path / "inference.treefile", sorted(sequences))
    assert descriptors["median_pair_distance"] == pytest.approx(other["median_pair_distance"])
