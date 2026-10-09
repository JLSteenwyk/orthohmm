"""Separate comparator routing preserves the official statistic and missing endpoints."""

import copy
import csv
import io
import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools import bootstrap_native_qfo_comparators as current
from benchmark_tools.audit_qfo_swiss_counts import REFERENCE_SHA
from benchmark_tools.bootstrap_qfo_swiss_comparators import METHODS


ROOT = Path(__file__).resolve().parents[2]


def literal_scores(counts):
    ppv = (counts["TP"] + 2) / (counts["TP"] + counts["FP"] + 4)
    tpr = (counts["TP"] + 2) / (counts["TP"] + counts["FN"] + 4)
    return dict(PPV=ppv, TPR=tpr, F1=2 * ppv * tpr / (ppv + tpr))


def fixture():
    families = [f"F{i:02d}" for i in range(18)]
    legacy = dict(status="corrected_comparison_swiss_counts_verified", families=families,
        reference=dict(sha256=REFERENCE_SHA), reference_relation_count=270,
        shared_represented_genes={}, publication_ready=False, uncertainty_admitted=False, methods=[])
    for index, method in enumerate(METHODS):
        records = []
        for i, family in enumerate(families):
            tp, fp = 1 + (i + index) % 7, (i + 2 * index) % 7
            counts = dict(TP=tp, FP=fp, FN=8-tp, TN=7-fp)
            records.append(dict(family=family, counts_without_prior=counts,
                represented_genes=[f"{family}_g{j}" for j in range(6)], statistics_with_prior=literal_scores(counts)))
        ppv = float(np.mean([r["statistics_with_prior"]["PPV"] for r in records]))
        tpr = float(np.mean([r["statistics_with_prior"]["TPR"] for r in records]))
        legacy["methods"].append(dict(method=method, status="counts_verified", families=records,
            aggregate=dict(PPV=ppv, TPR=tpr, F1=2*ppv*tpr/(ppv+tpr))))
    native = {cell: dict(copy.deepcopy(legacy["methods"][index]), cell=cell)
              for cell, index in zip(("p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p1_c0_r1"), (0, 1, 4, 6))}
    return native, legacy


def test_literal_paired_calculation_matches_all_metrics_and_48_endpoint_adjustment():
    native, legacy = fixture()
    before = copy.deepcopy((native, legacy))
    result = current.analyze(native, legacy, replicates=1000)
    assert (native, legacy) == before
    assert len(result["comparisons"]) == 16
    assert result["planned_endpoints"] == result["multiplicity_endpoints"] == 48
    assert result["estimated_endpoints"] == 24
    assert result["protocol_controls_match"] is False
    weights = np.random.Generator(np.random.PCG64(current.SEED)).multinomial(18, np.full(18, 1/18), size=1000)
    all_rows = {r["method"]: r for r in legacy["methods"]}
    all_rows.update(native)
    for comparison in result["comparisons"]:
        if comparison["candidate"] not in native:
            assert comparison["metrics"] is comparison["family_differences"] is None
            continue
        candidate, reference = (all_rows[comparison[key]] for key in ("candidate", "reference"))
        arrays = [np.array([[r["statistics_with_prior"][key] for key in ("PPV", "TPR")]
                           for r in item["families"]]) for item in (candidate, reference)]
        estimates, points, per_family = [], [], []
        for array in arrays:
            means = weights @ array / 18
            estimates.append(np.column_stack((2*means[:,0]*means[:,1]/means.sum(axis=1), means)))
            p, r = array.mean(axis=0)
            points.append(np.array([2*p*r/(p+r), p, r]))
            per_family.append(np.column_stack((2*array[:,0]*array[:,1]/array.sum(axis=1), array)))
        delta = estimates[0] - estimates[1]
        family_delta = per_family[0] - per_family[1]
        for j, metric in enumerate(current.METRICS):
            actual = comparison["metrics"][metric]
            assert actual["difference"] == pytest.approx(points[0][j]-points[1][j], abs=1e-14)
            assert actual["paired_percentile_ci"] == pytest.approx(np.quantile(delta[:,j], [.025,.975], method="linear"))
            assert actual["bonferroni_percentile_ci"] == pytest.approx(np.quantile(delta[:,j], [.05/96,1-.05/96], method="linear"))
            assert actual["family_wins"] == np.sum(family_delta[:,j] > 1e-10)
            assert actual["family_losses"] == np.sum(family_delta[:,j] < -1e-10)
            assert actual["family_wins"]+actual["family_ties"]+actual["family_losses"] == 18
    assert result["publication_ready"] is result["independent_confirmation"] is False


@pytest.mark.parametrize("change", ["empty", "cell", "key", "family", "gene", "truth", "negative",
    "boolean", "statistic", "aggregate", "nan", "legacy_scope", "replicates", "seed"])
def test_changed_counts_and_scope_refused(change):
    native, legacy = fixture()
    replicates, seed = 100, current.SEED
    row = native["p0_c0_r0"]
    if change == "empty": native = {}
    elif change == "cell": row["cell"] = "p1_c1_r1"
    elif change == "key": native["unknown"] = row
    elif change == "family": row["families"][0]["family"] = "missing"
    elif change == "gene": row["families"][0]["represented_genes"][0] = "unknown"
    elif change == "truth": row["families"][0]["counts_without_prior"]["FN"] += 1
    elif change == "negative": row["families"][0]["counts_without_prior"]["TP"] = -1
    elif change == "boolean": row["families"][0]["counts_without_prior"]["TP"] = True
    elif change == "statistic": row["families"][0]["statistics_with_prior"]["PPV"] += .01
    elif change == "aggregate": row["aggregate"]["F1"] += .01
    elif change == "nan": row["aggregate"]["F1"] = float("nan")
    elif change == "legacy_scope": legacy["publication_ready"] = True
    elif change == "replicates": replicates = 99
    elif change == "seed": seed = -1
    with pytest.raises(ValueError): current.analyze(native, legacy, replicates=replicates, seed=seed)


def test_missing_endpoint_table_never_imputes_or_shrinks_correction():
    native, legacy = fixture()
    result = current.analyze(native, legacy, replicates=100)
    tsv, markdown = current.tables(result)
    rows = list(csv.DictReader(io.StringIO(tsv), delimiter="\t"))
    assert len(rows) == 48
    assert sum(row["difference"] == "NA" for row in rows) == 24
    assert "48 planned endpoints" in markdown
    assert result["comparisons"][8]["missing_reason"].startswith("not_in_supplied_fresh_cell_snapshot")


def test_real_provenance_and_counts_prepare_without_production_bootstrap():
    native, legacy, provenance = current.prepare(ROOT)
    matrices = current.native_values(native, legacy)
    assert len(native) == 4 and len(matrices) == 6
    assert all(value.shape == (18, 2) for value in matrices.values())
    assert len(provenance["evidence"]) == 8
    assert provenance["protocol"]["sha256"] == current.PROTOCOL_SHA
    assert native["p0_c0_r1"]["resources"] is None


def test_run_refuses_existing_output_before_calculation(tmp_path, monkeypatch):
    output = tmp_path / "existing"
    output.mkdir()
    monkeypatch.setattr(current, "prepare", lambda root: pytest.fail("must refuse before preparing"))
    with pytest.raises(FileExistsError): current.run(ROOT, output)


def test_run_serializes_fixture_only_not_real_intervals(tmp_path, monkeypatch):
    native, legacy = fixture()
    result = current.analyze(native, legacy, replicates=100)
    result["protocol_controls_match"] = True
    monkeypatch.setattr(current, "prepare", lambda root: (native, legacy,
        dict(evidence=[], helpers=[], protocol=dict(path=str(Path(__file__)), **current.identity(Path(__file__))))))
    monkeypatch.setattr(current, "analyze", lambda *args: result)
    output = tmp_path / "fixture"
    observed = current.run(ROOT, output)
    assert observed["estimated_endpoints"] == 24
    assert json.loads((output / "report.json").read_text())["families"] == legacy["families"]
