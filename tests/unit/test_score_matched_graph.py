from copy import deepcopy
import json

import numpy as np
import pytest

from benchmark_tools.score_matched_graph import CONDITIONS, SEEDS, render, score_cell, summarize
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def rows():
    result = []
    for condition in CONDITIONS:
        for index, seed in enumerate(SEEDS):
            for arm in ("hmm", "diamond"):
                tp = 50 + index if arm == "hmm" else 50
                result.append(dict(condition=condition, seed=seed, arm=arm, truth={"same": True},
                    score=dict(tp=tp, fp=100-tp, fn=100-tp, input_genes=100, eligible_true_pairs=100,
                               f1=tp/100, precision=tp/100, recall=tp/100, undefined_ratios=[])))
    return result


def test_equal_seed_and_condition_means():
    result = summarize(rows())
    for condition in (*CONDITIONS, "overall"):
        r = result["contrasts"][condition]["f1"]
        assert r["difference_percentage_points"] == pytest.approx(2)
        assert r["hmm_mean"] == pytest.approx(.52)
        assert (r["wins"], r["ties"], r["losses"]) == (4, 1, 0)
        assert r["bonferroni_8_ci"][0] <= r["marginal_95_percent_ci"][0]
        assert r["bonferroni_8_ci"][1] >= r["marginal_95_percent_ci"][1]


def test_same_seed_blocks_carried_across_conditions():
    result = summarize(rows())
    assert result["contrasts"]["overall"]["f1"]["bonferroni_8_ci"] == pytest.approx(
        result["contrasts"]["baseline"]["f1"]["bonferroni_8_ci"])
    rng = np.random.Generator(np.random.PCG64(20260927))
    weights = rng.multinomial(5, np.full(5, .2), size=20000)
    expected = np.quantile(weights @ np.arange(5) / 5, [.025/8, 1-.025/8])
    assert result["contrasts"]["overall"]["f1"]["bonferroni_8_ci"] == pytest.approx(expected)


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "truth", "count", "calibration", "metric"])
def test_reject_nonprotocol_rows(mutation):
    values = rows()
    if mutation == "missing":
        values.pop()
    elif mutation == "duplicate":
        values.append(deepcopy(values[0]))
    elif mutation == "truth":
        values[0]["truth"] = {"different": True}
    elif mutation == "count":
        values[0]["score"]["input_genes"] += 1
    elif mutation == "calibration":
        values[0]["seed"] = 20261101
    else:
        values[0]["score"]["f1"] = .999
    with pytest.raises(ValueError):
        summarize(values)


def test_generated_table():
    result = summarize(rows())
    result["limitations"] = ["fixture"]
    text = render(result)
    assert "52.0000 | 50.0000 | +2.0000" in text
    assert "4/1/0" in text and "fixture" in text


def test_pairs_use_orthology_truth_not_ancestral_family(tmp_path):
    a, b = tmp_path / "s1.fasta", tmp_path / "s2.fasta"
    a.write_text(">a\nACD\n")
    b.write_text(">b\nACD\n>c\nACD\n")
    truth = tmp_path / "truth.json"
    truth.write_text(json.dumps(dict(extant_genes=3, ortholog_pairs=[["a", "b"], ["a", "c"]], families={"all": ["a", "b", "c"]})))
    prediction = tmp_path / "final.tsv"
    prediction.write_text("a\tb\nc\n")
    cell = dict(condition="baseline", seed=SEEDS[0], arm="hmm", genes=3, truth=record(truth),
                final_partition=record(prediction), execution=record(truth), native_receipt=record(truth), inputs=[record(a), record(b)])
    result = score_cell(cell)
    assert (result["score"]["tp"], result["score"]["fp"], result["score"]["fn"]) == (1, 0, 1)
    assert result["score"]["f1"] == pytest.approx(2/3)
    assert result["gene_coverage"] == 1 and result["nonsingleton_coverage"] == pytest.approx(2/3)
