import copy

import pytest

from benchmark_tools import prepare_ob_error_strata as strata


def family(size=2, identity=.5):
    return {"genes": size, "residue_lengths": [100] * size, "composition_flags": [False] * size,
            "species_counts": {f"s{i}": 1 for i in range(size)}, "mean_pairwise_identity": identity}


@pytest.mark.parametrize("size,expected", [(20, "small_2_20"), (21, "medium_21_50"), (50, "medium_21_50"), (51, "large_gt_50")])
def test_fixed_size_boundaries(size, expected):
    result = strata.assign_strata({"RefOG001.txt": family(size)})
    assert result["assignments"]["RefOG001.txt"]["size"] == expected
    assert len(result["strata"]) == 14 and result["multiplicity_endpoints"] == 84


def test_ties_missingness_and_order_independence():
    families = {"a": family(identity=.2), "b": family(identity=.5), "c": family(identity=.5),
                "d": family(identity=.9), "e": family(identity=None)}
    families["a"]["residue_lengths"] = [49, 151]
    families["b"]["residue_lengths"] = [50, 150]
    families["c"]["residue_lengths"] = [0, 100]
    families["a"]["composition_flags"] = [None, True]
    families["b"]["composition_flags"] = [None, False]
    families["a"]["species_counts"] = {"s": 2}
    result = strata.assign_strata(families)
    assert result == strata.assign_strata(dict(reversed(list(families.items()))))
    assert result["identity_median"] == .5
    assert result["strata"]["identity:lower_identity"] == ["a", "b", "c"]
    assert result["strata"]["identity:higher_identity"] == ["d"]
    assert result["strata"]["identity:missing"] == ["e"]
    assert result["assignments"]["a"]["relative_length"] == "short_relative"
    assert result["assignments"]["b"]["relative_length"] == "not_short_relative"
    assert result["assignments"]["c"]["relative_length"] == "missing"
    assert result["assignments"]["a"]["composition"] == "concentrated"
    assert result["assignments"]["b"]["composition"] == "missing"
    assert result["assignments"]["a"]["copy_number"] == "multi_copy"
    for dimension, categories in strata.CATEGORIES.items():
        members = [name for category in categories for name in result["strata"][dimension + ":" + category]]
        assert sorted(members) == sorted(families)


def test_all_missing_identities_keep_empty_observed_bins():
    result = strata.assign_strata({"a": family(identity=None)})
    assert result["identity_median"] is None
    assert result["strata"]["identity:lower_identity"] == []
    assert result["illustrative_families_by_stratum"]["identity:lower_identity"] is None


@pytest.mark.parametrize("problem", ["nonfinite", "flag", "length", "copy"])
def test_invalid_feature_inventory_rejected(problem):
    row = family()
    if problem == "nonfinite":
        row["mean_pairwise_identity"] = float("nan")
    elif problem == "flag":
        row["composition_flags"][0] = "False"
    elif problem == "length":
        row["residue_lengths"].pop()
    elif problem == "copy":
        row["species_counts"] = {"s": 3}
    with pytest.raises(ValueError):
        strata.assign_strata({"a": row})


@pytest.mark.parametrize("problem", [None, "failed", "missing", "source", "job"])
def test_whole_alignment_panel_required(problem):
    source = {"path": "/frozen/source.py"}
    panel = {"status": "reference_alignments_prepared_unscored", "accuracy_evaluated": False,
             "failed_families": [], "preflight": {"job_id": "21309", "status": "ready_unscored",
                 "accuracy_evaluated": False, "families": ["a", "b"], "source": copy.deepcopy(source)},
             "runs": [{"refog": name, "status": "alignment_validated", "exit_code": 0,
                       "accuracy_evaluated": False} for name in ("a", "b")]}
    if problem == "failed":
        panel["runs"][0]["status"] = "failed"
    elif problem == "missing":
        panel["runs"].pop()
    elif problem == "source":
        panel["preflight"]["source"] = {}
    elif problem == "job":
        panel["preflight"]["job_id"] = "wrong"
    if problem:
        with pytest.raises(ValueError):
            strata.check_alignment_panel(panel, {"a": set(), "b": set()}, source)
    else:
        strata.check_alignment_panel(panel, {"a": set(), "b": set()}, source)


def test_live_job_blocks_feature_and_reference_access(tmp_path, monkeypatch):
    monkeypatch.setattr(strata.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed\n21308|COMPLETED|0:0|00:00:23\n21309|RUNNING|0:0|00:01:00\n")
    monkeypatch.setattr(strata, "read_frozen", lambda *a: pytest.fail("Read features before job admission"))
    with pytest.raises(ValueError):
        strata.prepare(tmp_path, tmp_path / "out")
    assert not (tmp_path / "out").exists()
