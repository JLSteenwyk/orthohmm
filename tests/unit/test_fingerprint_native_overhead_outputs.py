import copy
import json
import itertools
from pathlib import Path

import pytest

from benchmark_tools import fingerprint_native_overhead_outputs as module

OWNERS = {"a": "sp1", "b": "sp2", "c": "sp1", "d": "sp3"}


def test_group_labels_and_order_are_not_biological_differences():
    left = {"OG0": ["b", "a"], "OG1": ["d", "c"]}
    right = {"renamed": ["c", "d"], "another": ["a", "b"]}
    assert module.partition_identity(left, OWNERS) == module.partition_identity(right, OWNERS)
    changed = {"OG0": ["d", "a"], "OG1": ["b", "c"]}
    assert module.partition_identity(left, OWNERS) != module.partition_identity(changed, OWNERS)


@pytest.mark.parametrize("groups", [{"a": ["a", "b"]}, {"a": ["a", "b", "c", "d", "e"]},
    {"a": ["a", "a", "b", "c", "d"]}, {"a": ["a", "b"], "b": ["b", "c", "d"]},
    {"a": [], "b": ["a", "b", "c", "d"]}])
def test_invalid_partitions_rejected(groups):
    with pytest.raises(ValueError):
        module.partition_identity(groups, OWNERS)


def test_pair_sets_ignore_order_orientation_and_redundant_rows():
    left = [("a", "b"), ("b", "c")]
    right = [("c", "b"), ("b", "a"), ("a", "b")]
    assert module.pair_identity(left, OWNERS) == module.pair_identity(right, OWNERS)
    assert module.pair_identity(left, OWNERS) != module.pair_identity([("a", "b"), ("a", "d")], OWNERS)
    assert module.pair_identity([], OWNERS)["pairs"] == 0


@pytest.mark.parametrize("pairs", [[("a", "a")], [("a", "c")], [("a", "unknown")],
                                    [("a",)], ["ab"], [("a", 1)]])
def test_invalid_pairs_rejected(pairs):
    with pytest.raises(ValueError):
        module.pair_identity(pairs, OWNERS)


def test_json_row_encoding_preserves_delimiter_and_unicode_boundaries():
    assert module.digest_rows([("a,b", "c")]) != module.digest_rows([("a", "b,c")])
    assert module.digest_rows([("a\nb", "c")]) != module.digest_rows([("a", "b"), ("c",)])
    assert len(module.digest_rows([("\u03b1", "\u03b2")])) == 64


@pytest.mark.parametrize("method", ["orthohmm_high_sensitivity", "orthohmm_satellite_v2"])
def test_native_orthohmm_fingerprint_uses_existing_strict_adapters(tmp_path, monkeypatch, method):
    out = tmp_path / "output"
    out.mkdir()
    (out / "orthohmm_orthogroups.txt").write_text("OG0: a b\nOG1: c d\n")
    phylo = out / "orthohmm_phylogeny"
    phylo.mkdir()
    (phylo / "orthohmm_root_hogs.tsv").write_text("root_hog\tsource_family\tgenes\nr1\tf1\tb,a\nr2\tf2\td,c\n")
    (phylo / "orthohmm_pairwise_orthologs.tsv").write_text("gene_a\tspecies_a\tgene_b\tspecies_b\na\tsp1\tb\tsp2\nc\tsp1\td\tsp3\n")
    run = {"native_method": method, "configuration": {"output": str(out)}, "dataset": {"inputs": []}}
    monkeypatch.setattr(module, "relocate_evidence", lambda run, roots: (copy.deepcopy(run), str))
    monkeypatch.setattr(module, "input_universe", lambda dataset: (OWNERS, ["sp1", "sp2", "sp3"]))
    result = module.fingerprint(run, {})
    assert result["identity"]["orthogroups"]["groups"] == 2
    assert result["scientific_timings_admitted"] is False
    assert len(result["evidence"]) == (3 if method.endswith("satellite_v2") else 1)
    if method.endswith("satellite_v2"):
        assert result["identity"]["root_hogs"] == result["identity"]["orthogroups"]
        assert result["identity"]["native_pairs"]["pairs"] == 2


@pytest.mark.parametrize("fault", [None, "missing_orientation", "disagreeing_orientation", "unknown_gene"])
def test_portable_orthofinder_checkpoint_and_native_pairs(tmp_path, monkeypatch, fault):
    out = tmp_path / "results"
    out.mkdir()
    (out / "SequenceIDs.txt").write_text("0_0: a\n0_1: c\n1_0: b\n2_0: d\n")
    (out / "clusters_OrthoFinder_I1.5.txt_id_pairs.txt").write_text("begin\n0 0_0 1_0 $\n1 0_1 2_0 $\n)\n")
    relations = {("sp1", "sp2"): ("a", "b"), ("sp2", "sp1"): ("b", "a"),
                 ("sp1", "sp3"): ("c", "d"), ("sp3", "sp1"): ("d", "c")}
    for a, b in itertools.permutations(["sp1", "sp2", "sp3"], 2):
        directory = out / "Orthologues" / ("Orthologues_" + a)
        directory.mkdir(parents=True, exist_ok=True)
        if fault == "missing_orientation" and (a, b) == ("sp2", "sp1"):
            continue
        pair = relations.get((a, b))
        if fault == "disagreeing_orientation" and (a, b) == ("sp2", "sp1"):
            pair = ("b", "c")
        if fault == "unknown_gene" and (a, b) == ("sp1", "sp2"):
            pair = ("unknown", "b")
        text = f"Orthogroup\t{a}\t{b}\n"
        if pair:
            text += f"OG0\t{pair[0]}\t{pair[1]}\n"
        (directory / f"{a}__v__{b}.tsv").write_text(text)
    run = {"native_method": "orthofinder_full", "configuration": {"output": str(out)}, "dataset": {"inputs": []}}
    monkeypatch.setattr(module, "relocate_evidence", lambda run, roots: (copy.deepcopy(run), str))
    monkeypatch.setattr(module, "input_universe", lambda dataset: (OWNERS, ["sp1", "sp2", "sp3"]))
    if fault:
        with pytest.raises(ValueError):
            module.fingerprint(run, {})
    else:
        result = module.fingerprint(run, {})
        assert result["identity"]["checkpoint_groups"] == module.partition_identity({"x": ["a", "b"], "y": ["c", "d"]}, OWNERS)
        assert result["identity"]["native_pairs"] == module.pair_identity([("a", "b"), ("c", "d")], OWNERS)
        assert len(result["evidence"]) == 8


def test_changed_group_file_during_fingerprinting_is_rejected(tmp_path, monkeypatch):
    path = tmp_path / "orthohmm_orthogroups.txt"
    path.write_text("OG0: a b\nOG1: c d\n")
    run = {"native_method": "orthohmm_high_sensitivity", "configuration": {"output": str(tmp_path)}, "dataset": {"inputs": []}}
    monkeypatch.setattr(module, "relocate_evidence", lambda run, roots: (copy.deepcopy(run), str))
    monkeypatch.setattr(module, "input_universe", lambda dataset: (OWNERS, ["sp1", "sp2", "sp3"]))
    original = module.read_predictions
    def changed(path, format):
        result = original(path, format)
        path.write_text("OG0: a b c d\n")
        return result
    monkeypatch.setattr(module, "read_predictions", changed)
    with pytest.raises(ValueError):
        module.fingerprint(run, {})


@pytest.mark.parametrize("index", [0, 1, 2])
def test_completed_small_fixture_native_outputs_can_be_fingerprinted(index):
    root = Path(__file__).resolve().parents[2]
    archive = root / "benchmarks/work/dgx_frontier_native_21831"
    if not archive.exists():
        pytest.skip("Retained raw engineering archive is not present")
    run = json.loads((archive / f"frontier_native_smoke_v1/run_{index:02d}/preparation.json").read_text())["run"]
    result = module.fingerprint(run, {Path("/home/jlsteenwyk/projects/orthohmm-publication"): archive})
    key = "checkpoint_groups" if index == 2 else "orthogroups"
    assert result["identity"][key]["genes"] == 645
    assert result["identity"][key]["groups"] == (99 if index == 2 else 98)
    if index:
        assert result["identity"]["native_pairs"]["pairs"] == (1834 if index == 2 else 1835)
