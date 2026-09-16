import pytest

from benchmark_tools.validate_factorial_native import check_native_metadata


def fixture(expanded):
    parameters = {"aligner": "mafft", "checkpoint_source": None, "cpu": 32,
                  "explicit_unconstrained_ablation": False, "pair_orthology_rule": "positive_paralogy",
                  "root_duplication_rule": "species_overlap", "species_tree": None,
                  "species_tree_mode": "infer", "species_tree_rooting": "min_variance",
                  "tree_builder": "FastTree",
                  "satellite_membership_policy": "high_confidence_pair" if expanded else "unconstrained"}
    cell = {"candidate_expansion": expanded, "argv": ["frozen", "command"]}
    metrics = {"status": "complete", "parameters": parameters, "command": cell["argv"][:]}
    native = {k: parameters[k] for k in ("pair_orthology_rule", "root_duplication_rule", "species_tree_mode", "species_tree_rooting")}
    native.update(mode="reconcile", cpu_budget=32,
                  species_tree_source="internally_inferred_from_orthohmm_single_copy_families",
                  membership_reconciliation={"policy": "high_confidence_pair", "constraints": 3,
                                             "supported_constraints": 2, "detached_constraints": 1} if expanded else None)
    return metrics, native, cell, 3 if expanded else 0


@pytest.mark.parametrize("expanded", [False, True])
def test_accepts_exact_native_policy(expanded):
    check_native_metadata(*fixture(expanded))


@pytest.mark.parametrize("change", ["command", "unfinished", "cpu", "supplied_tree", "rule", "policy", "constraints", "accounting", "checkpoint"])
def test_rejects_native_parameter_drift(change):
    metrics, native, cell, count = fixture(True)
    if change == "command":
        metrics["command"] = ["other"]
    elif change == "unfinished":
        metrics["status"] = "running"
    elif change == "cpu":
        native["cpu_budget"] = 64
    elif change == "supplied_tree":
        native["species_tree_source"] = "supplied"
    elif change == "rule":
        native["pair_orthology_rule"] = "other"
    elif change == "checkpoint":
        metrics["parameters"]["checkpoint_source"] = "/other"
    else:
        key = {"policy": "policy", "constraints": "constraints", "accounting": "detached_constraints"}[change]
        native["membership_reconciliation"][key] = "other" if key == "policy" else 10
    with pytest.raises(ValueError):
        check_native_metadata(metrics, native, cell, count)


def test_rejects_constraints_in_unexpanded_arm():
    metrics, native, cell, _ = fixture(False)
    with pytest.raises(ValueError):
        check_native_metadata(metrics, native, cell, 1)
