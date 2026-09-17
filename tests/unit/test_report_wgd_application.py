import pytest

from benchmark_tools.report_wgd_application import DIAGNOSTIC, METHODS, assemble


def fixture():
    owners = {"a": "Scerevisiae", "b": "Scerevisiae", "c": "Smikatae", "d": "Smikatae"}
    refs = {"p": list(owners)}
    row = {"orf_pair": ["a", "b"], "split_eligible": True, "reference_eligible": True,
           "reference_pillar": "p", "available_pillar_members": list(owners),
           "available_members_by_species": {"Scerevisiae": 2, "Smikatae": 2}, "experimental_class": "Sparse"}
    prepared = {"cohort_pairs": [row], "prespecified_examples": [{"orf_pair": ["a", "b"], "rank_sha256": "frozen"}]}
    admitted = {name: {"status": "admitted", "groups": {"x": ["a", "c"], "y": ["b", "d"]}}
                for name in (*METHODS, DIAGNOSTIC)}
    return prepared, refs, owners, admitted


def test_complete_rows_strata_examples_and_fixed_contrasts():
    report = assemble(*fixture())
    assert report["cohort_pairs"] == 1
    assert len(report["uncertainty"]["comparisons"]) == 12
    assert report["prespecified_examples"][0]["rank_sha256"] == "frozen"
    for entry in report["methods"].values():
        assert entry["summary"]["pairs"] == 1
        assert entry["strata"]["Sparse"]["pairs"] == 1
        assert entry["strata"]["Low"]["pairs"] == 0


def test_failures_remain_visible_without_zero_scores_or_replacement_examples():
    args = fixture()
    args[3]["sonicparanoid"] = {"status": "execution_failed", "reason": "timeout"}
    report = assemble(*args)
    assert report["methods"]["sonicparanoid"]["summary"] is None
    assert report["prespecified_examples"][0]["methods"]["sonicparanoid"] is None
    assert report["uncertainty"]["multiplicity"] == 12


@pytest.mark.parametrize("mutation", ["missing_method", "unfinished", "invented_failed_groups", "missing_example"])
def test_invalid_assembly_rejected(mutation):
    args = fixture()
    if mutation == "missing_method":
        del args[3][DIAGNOSTIC]
    elif mutation == "unfinished":
        args[3]["sonicparanoid"]["status"] = "running"
    elif mutation == "invented_failed_groups":
        args[3]["sonicparanoid"]["status"] = "execution_failed"
    else:
        args[0]["prespecified_examples"][0]["orf_pair"] = ["absent", "pair"]
    with pytest.raises(ValueError):
        assemble(*args)
