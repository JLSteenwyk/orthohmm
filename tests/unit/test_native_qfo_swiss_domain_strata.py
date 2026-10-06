import copy

import pytest

from benchmark_tools import export_native_qfo_swiss_domain_strata as export
from benchmark_tools import readback_native_qfo_swiss_domain_strata as reader


def panel():
    families = sorted(["APP", "BAR", "HOX", "NOX", "TRFE", "VATB", "MAPT", "PSEN"] +
                      [f"LOW{i}" for i in range(10)])
    memberships, genes = {}, {}
    for family in families:
        memberships[family] = [family + str(i) for i in range(4)]
        for i, gene in enumerate(memberships[family]):
            types = 2 if family in {"APP", "BAR", "HOX", "NOX", "TRFE", "VATB"} else 1
            domains = {f"pfam_{j}": {"instance": [[10 * j, 10 * j + 5, .01]]} for j in range(types)}
            if family in {"MAPT", "PSEN", "TRFE"} and i == 0:
                domains["pfam_0"]["instance"].append([25, 30, .02])
            if gene == "LOW00":
                domains = {}
            genes[gene] = dict(source_file="source.json", **export.annotations.features(dict(length=100, pfam=domains)))
    inventory = dict(genes=genes, families={f: export.annotations.summarize(m, genes) for f, m in memberships.items()})
    native, counts = dict(family_rows=[]), {cell: {} for cell in export.CELLS}
    for i, cell in enumerate(export.CELLS):
        for j, family in enumerate(families):
            c = dict(TP=20 + j - i, FP=(20 + 2 * j) * (1 - i), FN=5 + i * (j % 3), TN=10)
            counts[cell][family] = c
            native["family_rows"].append(dict(cell=cell, family=family, counts_without_prior=c,
                                             **export.sequence.family_statistics(c)))
    bins, rows, differences = export.project(native, memberships, inventory)
    report = dict(memberships=memberships, bins=bins, rows=rows, differences=differences,
                  family_rows=native["family_rows"])
    return native, memberships, inventory, report, counts


def test_complete_fixed_projection_independent_arithmetic():
    native, memberships, inventory, report, counts = panel()
    reader.verify_projection(report, counts, inventory)
    assert [len(v) for v in report["bins"].values()] == [18, 12, 6, 15, 3]
    assert len(report["rows"]) == 10
    assert len(report["differences"]) == 5
    assert len(report["family_rows"]) == 36
    assert report["differences"][0]["PPV"] > 0
    assert report["differences"][0]["TPR"] < 0
    means = sum(r["F1"] for r in native["family_rows"][:18]) / 18
    assert abs(means - report["rows"][0]["F1"]) > 1e-6
    assert "LOW00" in inventory["genes"]
    assert inventory["genes"]["LOW00"]["pfam_instance_count"] == 0


@pytest.mark.parametrize("mutation", ["missing_gene", "extra_gene", "shared_gene", "missing_family",
                                     "changed_summary", "wrong_bin", "missing_counts", "duplicate_counts",
                                     "wrong_counts", "negative_counts", "float_counts"])
def test_primary_refuses_changed_coverage_counts_or_bins(mutation):
    native, memberships, inventory, _, _ = panel()
    if mutation == "missing_gene":
        del inventory["genes"]["APP0"]
    elif mutation == "extra_gene":
        inventory["genes"]["extra"] = inventory["genes"]["APP0"]
    elif mutation == "shared_gene":
        memberships["APP"].append(memberships["BAR"][0])
    elif mutation == "missing_family":
        del memberships["APP"]
    elif mutation == "changed_summary":
        inventory["families"]["APP"]["annotated_genes"] -= 1
    elif mutation == "wrong_bin":
        for gene in memberships["APP"]:
            inventory["genes"][gene]["pfam_type_count"] = 0
        inventory["families"]["APP"] = export.annotations.summarize(memberships["APP"], inventory["genes"])
    elif mutation == "missing_counts":
        native["family_rows"].pop()
    elif mutation == "duplicate_counts":
        native["family_rows"].append(native["family_rows"][0])
    elif mutation == "wrong_counts":
        native["family_rows"][0]["F1"] = .1
    elif mutation == "negative_counts":
        native["family_rows"][0]["counts_without_prior"]["FN"] = -1
    else:
        native["family_rows"][0]["counts_without_prior"]["FN"] = 1.0
    with pytest.raises(ValueError):
        export.project(native, memberships, inventory)


@pytest.mark.parametrize("mutation", ["row_missing", "row_duplicate", "semantics", "member", "score",
                                     "difference", "bin", "family_count", "family_row_missing", "summary"])
def test_independent_reader_refuses_projection_changes(mutation):
    _, _, inventory, report, counts = panel()
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
        report["bins"][reader.NAMES[2]].pop()
    elif mutation == "family_count":
        report["family_rows"][0]["counts_without_prior"]["TN"] += 1
    elif mutation == "family_row_missing":
        report["family_rows"].pop()
    else:
        inventory["families"]["APP"]["annotated_genes"] -= 1
    with pytest.raises(ValueError):
        reader.verify_projection(report, counts, inventory)


@pytest.mark.parametrize("value", [
    dict(length=100, pfam={}),
    dict(length=100, pfam={"pfam_A": dict(instance=[[1, 20, .01], [30, 50, .01]])}),
    dict(length=100, pfam={"pfam_A": dict(instance=[[1, 20, .01]]),
                           "pfam_B": dict(instance=[[30, 50, .01]])}),
])
def test_independent_annotation_extraction(value):
    assert reader.annotation_features(value, "source.json") == dict(
        source_file="source.json", **export.annotations.features(value))


@pytest.mark.parametrize("value", [dict(length=True, pfam={}), dict(length=100),
    dict(length=100, pfam={"pfam_A": dict(instance=[[1, 101]])}),
    dict(length=100, pfam={"pfam_A": dict(instance=[[1, 20], [1, 20]])}),
    dict(length=100, pfam={"pfam_A": dict(instance=[])}),
    dict(length=100, pfam={"wrong": dict(instance=[[1, 20]])}),
])
def test_independent_annotation_refusals(value):
    with pytest.raises(ValueError):
        reader.annotation_features(value, "source.json")


def test_fresh_output_required_before_reading_sources(tmp_path):
    with pytest.raises(ValueError, match="Output already exists"):
        export.export(tmp_path / "missing_repo", tmp_path)


def test_dangling_output_symlink_refused(tmp_path):
    path = tmp_path / "result"
    path.symlink_to(tmp_path / "missing", target_is_directory=True)
    with pytest.raises(ValueError, match="Output already exists"):
        export.export(tmp_path / "missing_repo", path)
