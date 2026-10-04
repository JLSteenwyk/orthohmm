"""Complete-cohort and original-event guards for residual mechanism reporting."""

import dendropy
import pytest

from benchmark_tools import trace_simulation_oracle_residuals as trace
from orthohmm import phylogeny


def tree(text):
    return dendropy.Tree.get(data=text, schema="newick", preserve_underscores=True,
                             rooting="force-rooted")


def score(tp=1, fp=0, fn=0):
    return {"tp": tp, "fp": fp, "fn": fn, "predicted_pairs": tp + fp,
            "eligible_true_pairs": tp + fn, "f1": 2 * tp / (2 * tp + fp + fn),
            "precision": tp / (tp + fp), "recall": tp / (tp + fn)}


def report():
    return {"cells": [{"condition": c, "seed": s, "label": f"{c}_{s}",
            "status": "baseline_reproduced_oracle_scored", "candidates": [
                {"family": "Family0000000", "status": "oracle_eligible", "ancestral_families": ["1"],
                 "arms": {"generating_root": score(fn=int(s == 20261101))}},
                {"family": "Family0000001", "status": "unambiguous_bypass", "ancestral_families": ["2"],
                 "arms": {"generating_root": score(fp=int(s == 20261102))}}]}
            for c in sorted(trace.oracle.CONDITIONS) for s in range(20261101, 20261111)]}


def test_cohort_screens_all_cells_including_bypass_and_both_error_types():
    selected = trace.cohort(report())
    assert len(selected) == 14
    assert {c["status"] for _, c in selected} == {"oracle_eligible", "unambiguous_bypass"}


@pytest.mark.parametrize("change", ["missing", "duplicate_cell", "incomplete", "candidate_duplicate", "mixed", "unknown", "counts"])
def test_cohort_rejects_incomplete_or_changed_semantics(change):
    r = report()
    c = r["cells"][0]["candidates"][0]
    if change == "missing":
        r["cells"].pop()
    elif change == "duplicate_cell":
        r["cells"][-1] = r["cells"][0]
    elif change == "incomplete":
        r["cells"][0]["status"] = "failed"
    elif change == "candidate_duplicate":
        r["cells"][0]["candidates"].append(c)
    elif change == "mixed":
        c["ancestral_families"].append("9")
    elif change == "unknown":
        c["status"] = "mixed_ancestry_ineligible"
    else:
        c["arms"]["generating_root"]["fn"] += 1
    with pytest.raises(ValueError):
        trace.cohort(r)


def graph():
    return {"Root_1": ("D", ("Left_1", "Right_1")),
            "Left_1": ("S", ("a_1", "b_1")), "Right_1": ("S", ("c_1", "d_1")),
            **{g: ("F", ()) for g in ("a_1", "b_1", "c_1", "d_1")}}


def test_event_paths_reproduce_lca_not_leaf_order():
    leaves, pairs, paths, descendants = trace.history_index(graph())
    assert trace.common_ancestor(paths, "a_1", "b_1") == "Left_1"
    assert trace.common_ancestor(paths, "a_1", "d_1") == "Root_1"
    assert pairs == {("a_1", "b_1"), ("c_1", "d_1")}
    assert descendants["Root_1"] == leaves


@pytest.mark.parametrize("pair", [("a_1", "a_1"), ("a_1", "missing")])
def test_event_lca_requires_distinct_extant_genes(pair):
    with pytest.raises(ValueError, match="distinct"):
        trace.common_ancestor(trace.history_index(graph())[2], *pair)


def example(bypass=False, active=False, events=()):
    genes = {"F1__" + g for g in ("a_1", "b_1", "c_1", "d_1")}
    owners = {g: g.removeprefix("F1__").split("_")[0] for g in genes}
    truth = {("F1__a_1", "F1__b_1"), ("F1__c_1", "F1__d_1")}
    return trace.trace_candidate(phylogeny, tree("((a,b),(c,d));"),
        tree("((F1__a_1,F1__b_1),(F1__c_1,F1__d_1));"), graph(), "1", genes,
        owners, owners, truth, events, active, bypass, "Family0000000")


def test_true_duplication_with_reciprocal_losses_lacks_species_overlap():
    result = example()
    assert result["counts"] == {"tp": 2, "fp": 4, "fn": 0}
    assert result["error_classes"] == {"true_duplication_without_retained_species_overlap": 4}
    assert len(result["pairs"]) == 6
    errors = [r for r in result["pairs"] if r["error_class"]]
    assert all(r["history_event"] == "D" and r["pair_node"]["pair_event"] == "speciation" for r in errors)


def test_single_copy_bypass_keeps_predictions_not_new_oracle_pair_calls():
    result = example(bypass=True)
    assert result["counts"] == {"tp": 2, "fp": 4, "fn": 0}
    assert result["error_classes"] == {"single_copy_bypass_on_true_duplication": 4}


def test_constraint_supported_by_true_high_confidence_pair_is_not_detached():
    event = (7, {"source_genes": ["F1__a_1"], "target_genes": ["F1__b_1"]})
    result = example(active=True, events=[event])
    assert result["constraint_evidence"][0]["supported"] is True
    assert result["counts"]["fn"] == 0


@pytest.mark.parametrize("stage,expected", [("raw", "pair_rule_exclusion"),
    ("root", "root_partition_filter"), ("constraint", "unsupported_satellite_constraint")])
def test_false_negative_stage_precedence(stage, expected):
    row = {"truth": True, "predicted": False, "raw_predicted": stage != "raw",
           "same_root_group": stage != "root", "same_final_group": False}
    assert trace.error_class(row, False) == expected


def test_unexplained_or_inconsistent_false_negative_fails():
    with pytest.raises(ValueError, match="not explained"):
        trace.error_class({"truth": True, "predicted": False, "raw_predicted": True,
                           "same_root_group": True, "same_final_group": True}, False)


def test_parent_overlap_and_retained_overlap_distinguish_lost_signal():
    g = graph()
    g["Left_1"] = ("S", ("a_1", "b_1"))
    g["Right_1"] = ("S", ("a_2", "c_1"))
    del g["d_1"]
    g["a_2"] = ("F", ())
    owners = {"F1__a_1": "a", "F1__a_2": "a", "F1__b_1": "b", "F1__c_1": "c"}
    genes = set(owners) - {"F1__a_2"}
    result = trace.trace_candidate(phylogeny, tree("(a,(b,c));"), tree("((F1__a_1,F1__b_1),F1__c_1);"),
        g, "1", genes, {k: v for k, v in owners.items() if k in genes}, owners,
        {("F1__a_1", "F1__b_1")}, [], False, True, "Family0000000")
    r = next(r for r in result["pairs"] if r["gene_a"] == "F1__a_1" and r["gene_b"] == "F1__c_1")
    assert r["parent_species_overlap"] == ["a"]
    assert r["candidate_species_overlap"] == []


def test_new_event_reference_must_reproduce_retained_truth():
    genes = {"F1__" + g for g in ("a_1", "b_1", "c_1", "d_1")}
    owners = {g: g.removeprefix("F1__").split("_")[0] for g in genes}
    with pytest.raises(ValueError, match="local retained truth"):
        trace.trace_candidate(phylogeny, tree("((a,b),(c,d));"),
            tree("((F1__a_1,F1__b_1),(F1__c_1,F1__d_1));"), graph(), "1", genes,
            owners, owners, set(), [], False, False, "Family0000000")


def test_wrong_induced_topology_cannot_pass_merely_because_overlap_is_zero():
    genes = {"F1__" + g for g in ("a_1", "b_1", "c_1", "d_1")}
    owners = {g: g.removeprefix("F1__").split("_")[0] for g in genes}
    with pytest.raises(ValueError, match="original history clade"):
        trace.trace_candidate(phylogeny, tree("((a,b),(c,d));"),
            tree("((F1__a_1,F1__c_1),(F1__b_1,F1__d_1));"), graph(), "1", genes,
            owners, owners, {("F1__a_1", "F1__b_1"), ("F1__c_1", "F1__d_1")},
            [], False, False, "Family0000000")


def test_render_reports_all_candidates_and_no_independent_gain_claim():
    r = {"selection": "complete post hoc cohort", "screened_cells": 70, "screened_candidates": 10125,
         "candidates": [{"cell": "baseline_1", "family": "Family0000000", "status": "oracle_eligible",
                         "genes": ["a", "b"], "counts": {"tp": 1, "fp": 2, "fn": 3}, "error_classes": {"root_partition_filter": 3}}],
         "summary": {"error_classes": {"root_partition_filter": 3}}}
    rendered = trace.render(r)
    assert "10,125" in rendered and "baseline_1" in rendered and "| root_partition_filter | 3 |" in rendered
    assert "not independent accuracy gains" in rendered
