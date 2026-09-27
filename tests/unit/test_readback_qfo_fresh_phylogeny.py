import copy
import random

import pytest

from benchmark_tools.readback_qfo_fresh_phylogeny import (
    JOB, config_arguments, compare_pairs, pair_rows, verify_native, verify_summary,
)


def test_streaming_comparison_against_independent_sets():
    rng = random.Random(20260927)
    for _ in range(50):
        universe = [(f"a{i:02d}", f"z{i:02d}") for i in range(30)]
        a = {key: rng.choice(["high", "low"]) for key in universe if rng.random() < .6}
        b = {key: rng.choice(["high", "low"]) for key in universe if rng.random() < .6}
        result = compare_pairs(sorted(a.items()), sorted(b.items()))
        assert result["shared"] == len(a.keys() & b.keys())
        assert result["left_only"] == len(a.keys() - b.keys())
        assert result["right_only"] == len(b.keys() - a.keys())
        assert result["annotations_changed"] == sum(a[k] != b[k] for k in a.keys() & b.keys())
        assert result["pair_sets_equal"] == (a.keys() == b.keys())
        assert result["confidence_tables_equal"] == (a == b)


def test_empty_and_confidence_only_difference():
    assert compare_pairs([], [])["confidence_tables_equal"]
    result = compare_pairs([(("a", "b"), "high")], [(("a", "b"), "low")])
    assert result["pair_sets_equal"]
    assert not result["confidence_tables_equal"]


@pytest.mark.parametrize("data", ["b\tS1\ta\tS2\thigh\n", "a\tS1\tb\tS2\twrong\n",
    "a\tS1\tb\tS2\thigh\na\tS1\tb\tS2\thigh\n", "a\tS1\tb\n"])
def test_pair_reader_rejects_invalid_rows(tmp_path, data):
    path = tmp_path / "pairs"
    path.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\tconfidence\n" + data)
    with pytest.raises(ValueError):
        list(pair_rows(path))


def test_pair_reader_valid(tmp_path):
    path = tmp_path / "pairs"
    path.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\tconfidence\na\tS1\tb\tS2\thigh\n")
    assert list(pair_rows(path)) == [(("a", "b"), ("S1", "S2", "high"))]


def test_summary_schema_is_explicit():
    verify_summary(dict(schema_version=1, root_hogs=3), dict(root_hogs=3))
    for wrong in (dict(schema_version=2, root_hogs=3), dict(root_hogs=3),
                  dict(schema_version=1, root_hogs=4)):
        with pytest.raises(ValueError):
            verify_summary(wrong, dict(root_hogs=3))


def native_fixture():
    source = dict(path="/repo/benchmark_tools/run_qfo_fresh_phylogeny.py", bytes=1, sha256="source")
    plan = dict(repo="/repo", python="/python", environment={"PATH": "/bin"},
                aligner="/mafft", tree_builder="/FastTree", checked_records=[source])
    identity = dict(path="/plan", sha256="plan", bytes=2)
    start_record = dict(path="/started", sha256="started", bytes=3)
    started = dict(plan=identity, source=source, job_id=str(JOB), executable=plan["python"],
        environment=plan["environment"], config=config_arguments(plan), attempts=1, checkpoint_reuse=False)
    summary = dict(checkpoint_hits=0, remapped_checkpoint_hits=0,
                   species_tree_checkpoint_hit=False, reconciled_families=1, species_tree_families=2)
    complete = dict(status="fresh_phylogeny_complete_pending_readback", plan=identity,
                    started=start_record, accuracy_evaluated=False, summary=summary)
    return plan, identity, started, complete, start_record


def test_native_binding():
    verify_native(*native_fixture())


@pytest.mark.parametrize("index,key,value", [(2,"job_id","wrong"), (2,"executable","wrong"),
    (2,"environment",{}), (2,"config",{}), (2,"attempts",2), (2,"checkpoint_reuse",True),
    (2,"plan",{}), (3,"plan",{}), (3,"started",{}), (3,"accuracy_evaluated",True),
    (3,"status","incomplete")])
def test_native_mutations(index, key, value):
    args = copy.deepcopy(native_fixture())
    args[index][key] = value
    with pytest.raises(ValueError):
        verify_native(*args)


@pytest.mark.parametrize("key,value", [("checkpoint_hits",1), ("remapped_checkpoint_hits",1),
    ("species_tree_checkpoint_hit",True), ("reconciled_families",0), ("species_tree_families",0)])
def test_reject_reuse_or_missing_trees(key, value):
    args = native_fixture()
    args[3]["summary"][key] = value
    with pytest.raises(ValueError):
        verify_native(*args)
