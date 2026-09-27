from copy import deepcopy
import json

import numpy as np
import pytest

from benchmark_tools.trace_matched_graph_stages import (
    CONDITIONS, SEEDS, STAGES, aggregate, analyze_cell, component_pairs, edge_pairs, transition,
)
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.simulation_conditions import score_pairs


def test_direct_hits_collapse_directions_and_exclude_same_species():
    assert edge_pairs(["a", "b", "c"], dict(a=0, b=1, c=0),
                      [0, 1, 0, 0], [1, 0, 0, 2]) == {("a", "b")}


def test_components_keep_within_species_bridges_and_isolates():
    names, owners = ["a", "b", "c", "d"], dict(a=0, b=0, c=1, d=2)
    assert component_pairs(names, owners, [0, 1], [1, 2]) == {("a", "c"), ("b", "c")}
    assert component_pairs(names, owners, [], []) == set()


@pytest.mark.parametrize("sources,targets", [([0], []), ([-1], [0]), ([0], [4])])
def test_invalid_edges_rejected(sources, targets):
    with pytest.raises(ValueError):
        component_pairs(["a", "b"], dict(a=0, b=1), sources, targets)


def test_transitions_separate_true_false_and_gains_losses():
    truth = {("a", "b"), ("a", "c")}
    before = {("a", "b"), ("b", "c")}
    after = {("a", "c"), ("c", "d")}
    assert transition(before, after, truth) == dict(gained_true=1, lost_true=1, gained_false=1, lost_false=1)
    assert transition(before, before, truth) == dict(gained_true=0, lost_true=0, gained_false=0, lost_false=0)


def panel():
    return [dict(condition=c, seed=s, arm=a,
                 stages={stage: dict(tp=1, fp=1, fn=1, precision=.5, recall=.5, f1=.5) for stage in STAGES},
                 transitions={"initial_to_multipass": dict(gained_true=1)},
                 final_direct_support=dict(true_direct_supported=1))
            for c in CONDITIONS for s in SEEDS for a in ("hmm", "diamond")]


def test_equal_dataset_not_pooled_aggregation():
    rows = panel()
    rows[0]["stages"]["final"].update(tp=1000, f1=1.)
    result = aggregate(rows)
    assert result["baseline"]["hmm"]["stages"]["final"]["f1"] == .6
    assert result["overall"]["hmm"]["stages"]["final"]["f1"] == pytest.approx(.5 + .5/35)


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "wrong_seed", "wrong_arm"])
def test_panel_completeness(mutation):
    rows = deepcopy(panel())
    if mutation == "missing":
        rows.pop()
    elif mutation == "duplicate":
        rows.append(rows[0])
    elif mutation == "wrong_seed":
        rows[0]["seed"] = 1
    else:
        rows[0]["arm"] = "other"
    with pytest.raises(ValueError, match="complete paired"):
        aggregate(rows)


def test_full_cell_trace_and_frozen_score_gate(tmp_path):
    numeric = dict(gene_names=["a", "b", "c"], gene_to_species=[0, 1, 2],
                   hit_queries=[0, 1], hit_targets=[1, 2], hit_scores=[2., 3.])
    (tmp_path / "numeric.json").write_text(json.dumps(numeric))
    truth = [["a", "b"], ["a", "c"]]
    (tmp_path / "truth.json").write_text(json.dumps(dict(ortholog_pairs=truth)))
    outputs = []
    for stage in ("rbnh", "multipass"):
        path = tmp_path / (stage + "_edges.npz")
        np.savez(path, sources=np.array([0, 1]), targets=np.array([1, 2]), weights=np.array([2., 3.]))
        outputs.append(record(path))
    for stage in ("initial", "multipass", "final"):
        path = tmp_path / (stage + ".tsv")
        path.write_text("a\tb\nc\n" if stage == "initial" else "a\tb\tc\n")
        outputs.append(record(path))
    (tmp_path / "native.json").write_text(json.dumps(dict(outputs=outputs)))
    (tmp_path / "execution.json").write_text(json.dumps(dict(private_numeric=record(tmp_path / "numeric.json"))))
    cell = dict(execution=record(tmp_path / "execution.json"), native_receipt=record(tmp_path / "native.json"),
                final_partition=record(tmp_path / "final.tsv"), truth=record(tmp_path / "truth.json"),
                inputs=[], condition="baseline", seed=SEEDS[0], arm="hmm")
    expected = score_pairs([("a", "b"), ("a", "c"), ("b", "c")], truth, dict(a=0, b=1, c=2))
    result = analyze_cell(cell, expected)
    assert result["final_direct_support"] == dict(true_direct_supported=1, true_without_direct_hit=1,
                                                  false_direct_supported=1, false_without_direct_hit=0)
    assert result["transitions"]["initial_to_multipass"] == dict(gained_true=1, lost_true=0, gained_false=1, lost_false=0)
    assert result["stages"]["rbnh_components"]["tp"] == 2
    with pytest.raises(ValueError, match="frozen score"):
        analyze_cell(cell, dict(expected, tp=100))
    (tmp_path / "final.tsv").write_text("a\nb\nc\n")
    with pytest.raises(ValueError):
        analyze_cell(cell, expected)
