import math

import pytest

from benchmark_tools.prepare_corrected_swiss_sequence_strata import define_strata, descriptor_updates, CHANGED_SEQUENCES


def protein(entropy, length=100, canonical=None, fragment=False):
    return dict(length=length, canonical_residues=length if canonical is None else canonical,
                canonical_entropy_bits=None if entropy is None else entropy * math.log2(20),
                explicit_fragment_description=fragment)


def test_median_split_ties_and_all_secondary_bins():
    families = {"A": ["a"], "B": ["b"], "C": ["c"], "D": ["d"]}
    genes = {"a": protein(.7), "b": protein(.9), "c": protein(.9), "d": protein(1.)}
    result = define_strata(families, genes)
    assert result["primary_strata"] == {"lower_entropy": ["A", "B", "C"], "higher_entropy": ["D"], "missing_entropy": []}
    assert result["secondary_strata"]["concentrated"] == ["A"]
    assert result["secondary_strata"]["explicit_fragment"] == []


def test_no_silent_dropping_of_ineligible_member():
    result = define_strata({"F": ["a", "b"]}, {"a": protein(.9), "b": protein(.2, length=10)})
    assert result["primary_strata"]["missing_entropy"] == ["F"]
    assert result["secondary_strata"]["missing"] == ["F"]
    assert result["secondary_strata"]["short_relative"] == ["F"]


def test_true_flag_retained_despite_missing_other_member():
    result = define_strata({"F": ["a", "b"]}, {"a": protein(.7, fragment=True), "b": protein(None, canonical=0)})
    assert result["secondary_strata"]["concentrated"] == ["F"]
    assert result["secondary_strata"]["explicit_fragment"] == ["F"]
    assert result["median_family_entropy_cutoff"] is None


@pytest.mark.parametrize("row", [protein(float("nan")), protein(1.1), protein(None),
                                  protein(.9, canonical=101), protein(.9, length=0),
                                  protein(.9, canonical=89)])
def test_invalid_or_ineligible(row):
    if row["canonical_residues"] == 89:
        assert define_strata({"F": ["a"]}, {"a": row})["primary_strata"]["missing_entropy"] == ["F"]
    else:
        with pytest.raises(ValueError):
            define_strata({"F": ["a"]}, {"a": row})


def test_complete_reference_coverage_required():
    with pytest.raises(ValueError):
        define_strata({"F": ["a", "b"]}, {"a": protein(.9)})


@pytest.fixture
def updates():
    from copy import deepcopy
    names = [f"unchanged_{i}" for i in range(545)] + sorted(CHANGED_SEQUENCES) + ["A0A6I8Q293"]
    row = dict(canonical_entropy_bits=4., canonical_residues=100, length=100,
               maximum_canonical_frequency=.1, description="name PE=4 SV=2", source_file="Xenopus.fasta")
    old = {g: row.copy() for g in names}
    new = deepcopy(old)
    for gene in CHANGED_SEQUENCES:
        new[gene].update(canonical_entropy_bits=3.9, canonical_residues=90, length=90,
                         maximum_canonical_frequency=.2, description="name PE=4 SV=3")
    new["A0A6I8Q293"]["description"] = "name PE=3 SV=2"
    return old, new


def test_explicit_release_updates_only(updates):
    old, new = updates
    assert set(descriptor_updates(old, new)) == CHANGED_SEQUENCES | {"A0A6I8Q293"}


@pytest.mark.parametrize("fault", ["extra", "missing", "schema", "header", "fields"])
def test_unexpected_release_change_rejected(updates, fault):
    old, new = updates
    if fault == "extra":
        new["unchanged_0"]["length"] += 1
    elif fault == "missing":
        del new["unchanged_0"]
    elif fault == "schema":
        new["F6PXR2"]["extra"] = True
    elif fault == "header":
        new["A0A6I8Q293"]["description"] += " changed"
    else:
        new["F6PXR2"]["source_file"] = "other.fasta"
    with pytest.raises(ValueError):
        descriptor_updates(old, new)


def test_retained_corrected_inventory_replays_without_outcomes():
    import json
    from pathlib import Path
    from benchmark_tools.prepare_ob_candidate_neighborhood import check

    root = Path(__file__).resolve().parents[2]
    report = json.loads((root / "benchmark_tools/results/corrected_swiss_sequence_strata_20260918.json").read_text())
    assert report["prediction_statistics_evaluated"] is False
    assert report["publication_ready"] is False
    assert len(report["genes"]) == report["summary"]["matched_genes"] == 563
    assert report["summary"]["missing_genes"] == []
    assert len(report["recovered_accessions"]) == 14
    assert report["unchanged_original_descriptors"] == 545
    assert set(report["changed_original_descriptors"]) == CHANGED_SEQUENCES | {"A0A6I8Q293"}
    replay = define_strata(report["family_memberships"], report["genes"])
    assert replay == {k: report[k] for k in replay}
    assert [len(report["primary_strata"][k]) for k in ("lower_entropy", "higher_entropy", "missing_entropy")] == [9, 9, 0]
    assert report["secondary_strata"]["explicit_fragment"] == []
    assert report["secondary_strata"]["concentrated"] == []
    for key in ("protocol", "descriptor_update_protocol", "native_sequence_audit", "helper"):
        check(report[key])
