import pytest

from benchmark_tools.audit_unisave_fragment_source import choose_version
from benchmark_tools.verify_swiss_historical_fragments import selected_history, family_bins


def row(version, sv=1, first="12-Aug-2020", last="12-Aug-2020"):
    return dict(accession="P12345", sequenceVersion=sv, entryVersion=version,
                firstReleaseDate=first, lastReleaseDate=last)


@pytest.mark.parametrize("rows,sv", [([row(1)], 1),
    ([row(2, first="02-Dec-2020", last="02-Dec-2020"), row(1)], 1),
    ([row(3, 2, "10-Feb-2021", "10-Feb-2021"), row(2, 2, "02-Dec-2020", "02-Dec-2020"), row(1)], 2)])
def test_independent_selection_agrees_on_valid_histories(rows, sv):
    history = dict(results=rows)
    assert selected_history(history, "P12345", sv) == choose_version(history, "P12345", sv)


@pytest.mark.parametrize("rows", [[], [row(1), row(1)], [row(1), row(2)],
    [row(1, first="02-Dec-2020")], [row(1, first="17-Jun-2020", last="17-Jun-2020")],
    [dict(row(1), accession="OTHER")]])
def test_bad_history_rejected(rows):
    with pytest.raises(ValueError):
        selected_history(dict(results=rows), "P12345", 1)


def annotation(positive=False, feature=False, later=False):
    return dict(fragment_flag=positive, incomplete_sequence_features=[{}] if feature else [],
                selection_class="later_sequence_version" if later else "baseline_release")


def test_bins_preserve_missing_and_later_sensitivity():
    families = {"positive_missing": ["a", "b"], "unflagged": ["c"],
                "feature": ["d"], "later": ["e"], "missing": ["f", "c"]}
    annotations = dict(a=annotation(True), b=None, c=annotation(),
                       d=annotation(feature=True), e=annotation(True, later=True), f=None)
    bins = family_bins(families, annotations)
    assert bins == dict(annotation_positive=["feature", "later", "positive_missing"],
                        all_matched_unflagged=["unflagged"], missing_without_positive=["missing"])
    baseline = family_bins(families, annotations, True)
    assert baseline["annotation_positive"] == ["feature", "positive_missing"]
    assert baseline["missing_without_positive"] == ["later", "missing"]


@pytest.mark.parametrize("families,annotations", [({"f": []}, {}), ({"f": ["a", "a"]}, {"a": None}),
    ({"f": ["a"]}, {}), ({"f": ["a"]}, {"a": None, "extra": None})])
def test_incomplete_or_ambiguous_bins_rejected(families, annotations):
    with pytest.raises(ValueError):
        family_bins(families, annotations)
