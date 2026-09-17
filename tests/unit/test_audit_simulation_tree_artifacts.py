import pytest

from benchmark_tools.audit_simulation_tree_artifacts import index_panel, contrast, METHODS
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def panel():
    trees = [{"condition": str(c), "seed": s, "variant": v}
             for c in range(7) for s in range(10) for v in ("generating", "nni1", "nni2")]
    return {"status": "tree_panel_verified_unscored", "accuracy_evaluated": False,
            "records": [{**t, "method": m, "status": "admitted"} for t in trees for m in METHODS]}, trees


@pytest.mark.parametrize("problem", [None, "missing", "duplicate", "unknown", "scored", "unadmitted"])
def test_complete_inventory_gate(problem):
    report, trees = panel()
    if problem == "missing":
        report["records"].pop()
    elif problem == "duplicate":
        report["records"][-1] = report["records"][0]
    elif problem == "unknown":
        report["records"][0]["status"] = "running"
    elif problem == "scored":
        report["accuracy_evaluated"] = True
    elif problem == "unadmitted":
        report["status"] = "validating"
    if problem:
        with pytest.raises(ValueError):
            index_panel(report, trees)
    else:
        assert len(index_panel(report, trees)) == 420


@pytest.mark.parametrize("status", ["failed", "unavailable", "supplied_tree_not_retained"])
def test_unavailable_arms_never_equivalent(status):
    result = contrast(METHODS[0], {"status": status}, {"status": "admitted"})
    assert result["status"] == "unavailable"
    assert result["arm_statuses"]["target"] == status


def test_empty_inventory_rejected():
    arm = {"status": "admitted", "retained_artifacts": {}}
    with pytest.raises(ValueError, match="Empty"):
        contrast(METHODS[0], arm, arm)


@pytest.mark.parametrize("changed", [False, True])
def test_changed_gene_tree_is_retained_not_filtered(tmp_path, changed):
    a, b = tmp_path / "a", tmp_path / "b"
    a.write_text("(a,b);\n")
    b.write_text("(a,c);\n" if changed else "(a,b);\n")
    arms = [{"status": "admitted", "retained_artifacts": {"gene_tree": record(p)}} for p in (a, b)]
    result = contrast(METHODS[0], *arms)
    assert result["status"] == ("retained_upstream_different" if changed else "retained_upstream_equivalent")
    assert result["representation_exception_used"] is False


def test_mcl_command_path_exception_not_membership_exception(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    a.write_text("# cline: first\n0 1 2 $\n")
    b.write_text("# cline: second\n0 1 2 $\n")
    key = "clusters_OrthoFinder_I1.2.txt"
    arms = [{"status": "admitted", "retained_artifacts": {key: record(p)}} for p in (a, b)]
    result = contrast("orthofinder_full", *arms)
    assert result["status"] == "retained_upstream_equivalent"
    assert result["representation_exception_used"]
    b.write_text("# cline: second\n0 1 3 $\n")
    with pytest.raises(ValueError):
        contrast("orthofinder_full", *arms)
    arms[1]["retained_artifacts"][key] = record(b)
    assert contrast("orthofinder_full", *arms)["status"] == "retained_upstream_different"
