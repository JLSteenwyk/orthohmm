from itertools import combinations

import numpy as np
import pytest

from benchmark_tools.score_ygob_groups import paired_bootstrap, read_predictions, score_groups


def pairs(groups, universe):
    return {frozenset(pair) for group in groups.values()
            for pair in combinations(set(group) & universe, 2)}


def test_explicit_pair_enumeration_random_partitions():
    rng = np.random.default_rng(912)
    genes = [f"speciesA_{i}" for i in range(24)]
    references = {f"P{i}": genes[i * 4:(i + 1) * 4] for i in range(5)}
    references["singleton"] = [genes[20]]
    universe = set().union(*map(set, references.values()))
    truth = pairs(references, universe)
    for _ in range(100):
        predictions = {}
        for gene in genes:
            assignment = int(rng.integers(-1, 8))
            if assignment >= 0:
                predictions.setdefault(str(assignment), []).append(gene)
        score = score_groups(predictions, references, genes)
        predicted = pairs(predictions, universe)
        assert score["counts"] == {"tp": len(truth & predicted), "fp": len(predicted - truth),
                                   "fn": len(truth - predicted)}
        owners = {g: p for p, group in references.items() for g in group}
        for record in score["records"]:
            expected_fp = sum(sum(owners[g] == record["pillar"] for g in pair) / 2
                              for pair in predicted - truth)
            assert record["fp"] == expected_fp


def test_singletons_missing_genes_and_projection():
    refs = {"P": ["a", "b"], "S": ["c"], "T": ["d"]}
    score = score_groups({"one": ["a", "c", "excluded"], "two": ["b"]},
                         refs, {"a", "b", "c", "d", "excluded"})
    assert score["counts"] == {"tp": 0, "fp": 1, "fn": 1}
    assert [r["fp"] for r in score["records"]] == [0.5, 0.5, 0]
    assert score["reference_gene_coverage"] == 0.75
    assert score["excluded_predicted_genes"] == 1
    assert score["exact_reference_groups"] == 0
    assert score_groups({}, refs, "abcd")["counts"] == {"tp": 0, "fp": 0, "fn": 1}
    perfect = score_groups({"one": ["a", "b", "excluded"], "two": ["c"], "three": ["d"]},
                           refs, {*"abcd", "excluded"})
    assert perfect["metrics"] == {"f1": 1, "precision": 1, "recall": 1}
    assert perfect["exact_reference_groups"] == 3


@pytest.mark.parametrize("predictions", [{"P": ["a", "a"]}, {"P": ["a"], "Q": ["a"]},
                                        {"P": ["unknown"]}, {"P": []}])
def test_invalid_predictions(predictions):
    with pytest.raises(ValueError):
        score_groups(predictions, {"P": ["a", "b"]}, "ab")


def test_reference_validation():
    for refs in ({}, {"P": ["a", "a"]}, {"P": ["missing"]}):
        with pytest.raises(ValueError):
            score_groups({}, refs, "ab")


def test_prediction_adapters(tmp_path):
    path = tmp_path / "groups.txt"
    path.write_text("OG0: a b\nOG1: c\n")
    assert read_predictions(path, "named_groups") == {"OG0": ["a", "b"], "OG1": ["c"]}
    path.write_text("root_hog\tsource_family\tgenes\nR0\tF0\ta,b\n")
    assert read_predictions(path, "root_hogs") == {"R0": ["a", "b"]}
    path.write_text("OG0: a\nOG0: b\n")
    with pytest.raises(ValueError):
        read_predictions(path, "named_groups")
    with pytest.raises(ValueError):
        read_predictions(path, "root_hogs")


def test_bootstrap_matches_direct_resampling_and_batch_invariance():
    refs = {"P": ["a", "b", "c"], "Q": ["d", "e"], "S": ["f"]}
    base = score_groups({"A": ["a", "b", "d"], "B": ["c", "e", "f"]}, refs, "abcdef")
    other = score_groups(refs, refs, "abcdef")
    scores = {"base": base, "other": other, "identical": base}
    result = paired_bootstrap(scores, "base", replicates=200, seed=12, batch_size=17)
    assert result == paired_bootstrap(scores, "base", replicates=200, seed=12, batch_size=200)
    assert result["multiplicity_count"] == 6
    weights = np.random.Generator(np.random.PCG64(12)).multinomial(3, [1 / 3] * 3, size=200)
    differences = []
    for draw in weights:
        metrics = []
        for score in (other, base):
            tp, fp, fn = [sum(w * r[k] for w, r in zip(draw, score["records"]))
                          for k in ("tp", "fp", "fn")]
            metrics.append(2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0)
        differences.append(100 * (metrics[0] - metrics[1]))
    assert result["comparisons"]["other"]["f1"]["paired_95_percent_ci"] == pytest.approx(
        np.quantile(differences, [0.025, 0.975]))
    assert result["comparisons"]["identical"]["f1"]["paired_95_percent_ci"] == [0, 0]


def test_bootstrap_rejects_mismatched_references():
    one = score_groups({}, {"P": ["a"]}, "ab")
    two = score_groups({}, {"P": ["a", "b"]}, "ab")
    with pytest.raises(ValueError, match="share reference"):
        paired_bootstrap({"one": one, "two": two}, "one", replicates=100)
