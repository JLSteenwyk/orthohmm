import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools.analyze_swiss_domain_strata import CONTRASTS, METRICS, PRIMARY, SECONDARY, analyze, strata


def inputs():
    root = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    return (json.loads((root / "qfo_swiss_comparator_counts_20260917.json").read_text()),
            json.loads((root / "swiss_domain_annotation_inventory_20260917.json").read_text()))


def test_frozen_bins_and_descriptive_only_repeat_subset():
    counts, annotations = inputs()
    report = analyze(counts, annotations, replicates=100, seed=42)
    assert [len(report["strata"][n]) for n in (*PRIMARY, *SECONDARY)] == [12, 6, 15, 3]
    assert report["multiplicity_endpoints"] == 27
    for contrast in report["contrasts"]:
        assert len(contrast["family_differences"]) == 18
        assert all(set(m) == set(METRICS) for m in contrast["secondary_descriptive_differences"].values())


def test_direct_raw_count_resampling_reproduces_strata_and_interactions():
    counts, annotations = inputs()
    report = analyze(counts, annotations, replicates=100, seed=42)
    rng = np.random.Generator(np.random.PCG64(42))
    draws = {}
    for name in PRIMARY:
        families = report["strata"][name]
        weights = rng.multinomial(len(families), np.full(len(families), 1/len(families)), size=100)
        all_methods = []
        for method in counts["methods"]:
            by_family = {r["family"]: r["counts_without_prior"] for r in method["families"]}
            results = []
            for w in weights:
                selected = [by_family[f] for f, n in zip(families, w) for _ in range(n)]
                p = sum((r["TP"]+2)/(r["TP"]+r["FP"]+4) for r in selected)/len(selected)
                r = sum((r["TP"]+2)/(r["TP"]+r["FN"]+4) for r in selected)/len(selected)
                results.append([2*p*r/(p+r), p, r])
            all_methods.append(results)
        draws[name] = np.asarray(all_methods)
    for contrast, (candidate, reference) in zip(report["contrasts"], CONTRASTS):
        differences = {n: draws[n][candidate]-draws[n][reference] for n in PRIMARY}
        rows = [(r["metrics"], differences[r["stratum"]]) for r in contrast["primary_strata"]]
        rows.append((contrast["interaction_higher_minus_lower"], differences[PRIMARY[1]]-differences[PRIMARY[0]]))
        for metrics, delta in rows:
            for j, metric in enumerate(METRICS):
                assert metrics[metric]["nominal95"] == pytest.approx(np.quantile(delta[:, j], [.025, .975]))
                assert metrics[metric]["bonferroni27"] == pytest.approx(np.quantile(delta[:, j], [.05/54, 1-.05/54]))


@pytest.mark.parametrize("mutation", ["coverage", "summary", "family"])
def test_changed_annotation_inventory_rejected(mutation):
    counts, annotations = inputs()
    if mutation == "coverage":
        annotations["genes"].pop(next(iter(annotations["genes"])))
    elif mutation == "summary":
        annotations["families"]["APP"]["median_pfam_types_among_annotated"] = 1
    else:
        annotations["families"].pop("APP")
    with pytest.raises(ValueError):
        strata(counts, annotations)
