import copy
from fractions import Fraction

import pytest

from benchmark_tools import export_native_qfo_swiss_duplication_strata as export
from benchmark_tools import readback_native_qfo_swiss_duplication_strata as reader


def feature_panel(ratios=("0", "1/3", "2/3", "1", None)):
    memberships, families = {}, {}
    for i, ratio in enumerate(ratios):
        family = f"F{i}"
        memberships[family] = [family + str(j) for j in range(4)]
        value = Fraction(ratio) if ratio is not None else None
        numerator, denominator = (value.numerator, value.denominator) if value is not None else (0, 0)
        families[family] = dict(explicit_duplication_nodes=numerator, explicit_speciation_nodes=0,
                                default_speciation_nodes=denominator - numerator, informative_nodes=denominator,
                                child_overlap_nodes=0, mapped_members=4, duplication_fraction=ratio)
    values = sorted(Fraction(v) for v in ratios if v is not None)
    midpoint = (values[(len(values) - 1) // 2] + values[len(values) // 2]) / 2 if values else None
    groups = {name: [] for name in export.NAMES[1:]}
    for family, value in families.items():
        ratio = value["duplication_fraction"]
        groups[export.NAMES[3] if ratio is None else export.NAMES[1 + int(Fraction(ratio) > midpoint)]].append(family)
    return memberships, dict(families=families, median_fraction=str(midpoint) if midpoint is not None else None,
                             primary_strata=groups)


def panel():
    memberships, feature = feature_panel()
    counts, native = {cell: {} for cell in export.CELLS}, dict(family_rows=[])
    for i, cell in enumerate(export.CELLS):
        for j, family in enumerate(memberships):
            value = dict(TP=20 + j - i, FP=(20 + 9 * j) * (1 - i), FN=5 + i * (j % 3), TN=10)
            counts[cell][family] = value
            native["family_rows"].append(dict(cell=cell, family=family, counts_without_prior=value,
                                             **export.sequence.family_statistics(value)))
    bins, rows, differences = export.project(native, memberships, feature)
    report = dict(memberships=memberships, bins=bins, rows=rows, differences=differences,
                  family_rows=native["family_rows"])
    return native, memberships, feature, report, counts


def test_complete_projection_matches_independent_arithmetic():
    native, memberships, feature, report, counts = panel()
    reader.verify_projection(report, counts, feature)
    assert [len(v) for v in report["bins"].values()] == [5, 2, 2, 1]
    assert len(report["rows"]) == 8 and len(report["differences"]) == 4
    assert report["differences"][0]["PPV"] > 0 and report["differences"][0]["TPR"] < 0
    mean_family_f1 = sum(r["F1"] for r in native["family_rows"][:5]) / 5
    assert abs(mean_family_f1 - report["rows"][0]["F1"]) > 1e-6


@pytest.mark.parametrize("ratios,expected", [
    (("0", "1/3", "2/3", "1"), [4, 2, 2, 0]),
    (("0", "1/2", "1/2", "1"), [4, 3, 1, 0]),
    (("0", "1/2", "1"), [3, 2, 1, 0]),
    ((None, None), [2, 0, 0, 2]),
    (("1/2", "1/2", None), [3, 2, 0, 1]),
])
def test_exact_median_ties_and_missingness(ratios, expected):
    memberships, feature = feature_panel(ratios)
    primary = export.bins_for(memberships, feature)
    assert primary == reader.independent_bins(memberships, feature)
    assert [len(v) for v in primary.values()] == expected


def test_empty_bin_is_null_not_zero():
    native, memberships, feature, _, counts = panel()
    memberships.pop("F4")
    feature["families"].pop("F4")
    feature["primary_strata"][export.NAMES[3]] = []
    native["family_rows"] = [r for r in native["family_rows"] if r["family"] != "F4"]
    bins, rows, differences = export.project(native, memberships, feature)
    report = dict(memberships=memberships, bins=bins, rows=rows, differences=differences,
                  family_rows=native["family_rows"])
    reader.verify_projection(report, counts, feature)
    assert all(row[k] is None for row in [rows[3], rows[7], differences[3]] for k in export.METRICS)
    rows[3]["F1"] = 0
    with pytest.raises(ValueError):
        reader.verify_projection(report, counts, feature)


def test_alias_sensitive_informative_denominator_not_unique_genes_minus_one():
    memberships, feature = feature_panel(("3/8",))
    row = feature["families"]["F0"]
    row.update(explicit_duplication_nodes=21, default_speciation_nodes=35, informative_nodes=56,
               child_overlap_nodes=1, mapped_members=56)
    memberships["F0"] = [f"G{i}" for i in range(56)]
    assert export.bins_for(memberships, feature) == reader.independent_bins(memberships, feature)


@pytest.mark.parametrize("mutation", ["negative", "boolean", "float", "sum", "overlap", "mapped",
                                     "fraction", "noncanonical", "median", "tie_split", "missing_family",
                                     "shared_gene", "empty_genes", "duplicate_gene"])
def test_both_bin_readers_refuse_invalid_features(mutation):
    memberships, feature = feature_panel()
    row = feature["families"]["F1"]
    if mutation in ("negative", "boolean", "float"):
        row["explicit_duplication_nodes"] = {"negative": -1, "boolean": True, "float": 1.0}[mutation]
    elif mutation == "sum":
        row["default_speciation_nodes"] += 1
    elif mutation == "overlap":
        row["child_overlap_nodes"] = 4
    elif mutation == "mapped":
        row["mapped_members"] = 3
    elif mutation == "fraction":
        row["duplication_fraction"] = "1/2"
    elif mutation == "noncanonical":
        row["duplication_fraction"] = "2/6"
    elif mutation == "median":
        feature["median_fraction"] = "1/3"
    elif mutation == "tie_split":
        feature["primary_strata"][export.NAMES[1]].remove("F1")
        feature["primary_strata"][export.NAMES[2]].append("F1")
    elif mutation == "missing_family":
        del feature["families"]["F1"]
    elif mutation == "shared_gene":
        memberships["F1"][0] = memberships["F0"][0]
    elif mutation == "empty_genes":
        memberships["F1"] = []
    else:
        memberships["F1"][0] = memberships["F1"][1]
    for function in (export.bins_for, reader.independent_bins):
        with pytest.raises(ValueError):
            function(memberships, feature)


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "score", "negative", "boolean", "float"])
def test_primary_refuses_bad_native_count_vectors(mutation):
    native, memberships, feature, _, _ = panel()
    if mutation == "missing":
        native["family_rows"].pop()
    elif mutation == "duplicate":
        native["family_rows"].append(native["family_rows"][0])
    elif mutation == "score":
        native["family_rows"][0]["F1"] = .3
    else:
        native["family_rows"][0]["counts_without_prior"]["FN"] = {
            "negative": -1, "boolean": True, "float": 1.0}[mutation]
    with pytest.raises(ValueError):
        export.project(native, memberships, feature)


@pytest.mark.parametrize("mutation", ["row_missing", "row_duplicate", "semantics", "member", "score",
                                     "difference", "bin", "family_count", "family_row_missing", "status"])
def test_independent_reader_refuses_changed_results(mutation):
    _, _, feature, report, counts = panel()
    report = copy.deepcopy(report)
    if mutation == "row_missing":
        report["rows"].pop()
    elif mutation == "row_duplicate":
        report["rows"].append(report["rows"][0])
    elif mutation == "semantics":
        report["rows"][0]["prediction_semantics"] = "resolved_native_pairs"
    elif mutation == "member":
        report["rows"][0]["family_members"] = []
    elif mutation == "score":
        report["rows"][0]["F1"] += .1
    elif mutation == "difference":
        report["differences"][0]["TPR"] *= -1
    elif mutation == "bin":
        report["bins"][export.NAMES[2]].pop()
    elif mutation == "family_count":
        report["family_rows"][0]["counts_without_prior"]["TN"] += 1
    elif mutation == "family_row_missing":
        report["family_rows"].pop()
    else:
        report["rows"][0]["status"] = "empty_bin"
    with pytest.raises(ValueError):
        reader.verify_projection(report, counts, feature)


def test_existing_output_refused_before_any_sources(tmp_path):
    with pytest.raises(ValueError, match="Output already exists"):
        export.export(tmp_path / "missing", tmp_path)


def test_dangling_output_symlink_refused(tmp_path):
    output = tmp_path / "output"
    output.symlink_to(tmp_path / "missing", target_is_directory=True)
    with pytest.raises(ValueError, match="Output already exists"):
        export.export(tmp_path / "missing", output)


def mapped_panel():
    return (dict(HOX=["U1", "U2"]),
            dict(families=dict(HOX=dict(mapped_labels={"ENSEMBL_A": 10, "ENSEMBL_B": 20, "ENSEMBL_ALIAS": 20},
                                        exact_match=True, mapped_members=2))),
            dict(U1=10, U2=20))


def test_entry_join_preserves_alias_collision_and_distinct_accessions():
    for function in (export.check_memberships, reader.check_memberships):
        function(*mapped_panel())


@pytest.mark.parametrize("mutation", ["unknown", "wrong", "boolean", "negative", "native_collision",
                                     "missing_family", "retained_invalid", "count", "status"])
def test_entry_join_refuses_mismatch_not_just_equal_member_counts(mutation):
    memberships, mapping, identifiers = mapped_panel()
    if mutation == "unknown":
        del identifiers["U1"]
    elif mutation in ("wrong", "boolean", "negative", "native_collision"):
        identifiers["U1"] = {"wrong": 30, "boolean": True, "negative": -1, "native_collision": 20}[mutation]
    elif mutation == "missing_family":
        mapping["families"] = {}
    elif mutation == "retained_invalid":
        mapping["families"]["HOX"]["mapped_labels"]["ENSEMBL_A"] = True
    elif mutation == "count":
        mapping["families"]["HOX"]["mapped_members"] = 3
    else:
        mapping["families"]["HOX"]["exact_match"] = False
    for function in (export.check_memberships, reader.check_memberships):
        with pytest.raises(ValueError):
            function(memberships, mapping, identifiers)
