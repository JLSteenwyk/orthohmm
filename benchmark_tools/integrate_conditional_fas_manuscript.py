"""Integrate retained conditional FAS results without rewriting prior evidence."""

from itertools import combinations
import json
import math
from pathlib import Path

from benchmark_tools.audit_fas_stratum_weights import METHODS
from benchmark_tools.integrate_publication_supplements import record, require


PINS = {
    "parent": ("publication_supplements_manuscript_20261010_v1.md", "43daded442eab987cfe29863f8c661d82cdfbc4d55443d4e7896300e2917e45b"),
    "report": ("native_fas_conditional_report_20261010_v1/report.json", "83548e93afd5ffe0eeeef13250684dd65485b9efba1934d3377787659f7b3d7f"),
    "execution": ("native_fas_conditional_execution_20261010_v1.json", "fd3572b7a29f8cfbb7a9221ea5d92f6956bd6e9d344dff7eef87a90473d26675"),
}
ANCHORS = ("### Retrospective Native-Cell Comparator Sensitivity\n",
           "The completed simulations and tree perturbations do not cover arbitrary\n")


def sections(report, execution):
    require(report["status"] == "conditional_expected_native_ranges_complete" and
            report["all_eight_ranges_computed"] is True and report["conditional_design_ranges_computed"] is True,
            "Require all eight conditional ranges")
    require(report["alpha"] == .05 and report["component_error"] == .05/16 and
            [r["method"] for r in report["methods"]] == list(METHODS), "Changed allocation or methods")
    require(all(report[key] is False for key in ("observed_native_scores_changed","historical_scores_rerun",
            "unconditional_historical_interval_admission","biological_generalization_intervals",
            "other_endpoint_uncertainty_admitted","publication_ready")), "Changed conditional scope")
    require(execution["conditional_design_ranges_computed"] is True and all(execution[key] is False
            for key in ("unconditional_historical_interval_admission","biological_generalization_intervals","publication_ready")),
            "Changed execution scope")
    require(all(execution[key]["observations"][-1]["exit_code"] == 0
                for key in ("execution","tests","independent_recurrence_readback")) and
            execution["document_checks"][-1]["result"]["exit_code"] == 0, "Incomplete actual verification")
    require(len(report["contrasts"]) == 28, "Missing contrasts")
    for row in report["methods"]:
        require(row["status"] == "conditional_design_range_computed" and row["error"] is None and
                row["database_historically_hash_bound"] is False, "Incomplete range or changed historical scope")
        values = [row["observed_native_mean"]] + row["interval"]["expected_native_mean_bounds"]
        require(all(math.isfinite(x) and 0 <= x <= 1 for x in values) and values[1] <= values[2], "Invalid table values")
    for row,(left,right) in zip(report["contrasts"],combinations(report["methods"],2)):
        require((row["left"],row["right"]) == (left["method"],right["method"]), "Changed contrast identities")
        a,b = left["interval"]["expected_native_mean_bounds"],right["interval"]["expected_native_mean_bounds"]
        expected = [a[0]-b[1],a[1]-b[0]]
        require(row["conditional_expected_difference_bounds"] == expected and
                row["observed_difference"] == left["observed_native_mean"]-right["observed_native_mean"], "Changed projections")
    lines = ["### Conditional Sampling In Retained Eight-Method FAS", "",
        "This analysis concerns the retained September eight-method comparator snapshot,",
        "not the later four-cell fresh-native ablation above. It adds algorithmic",
        "sampling uncertainty without changing an observed score, method or default.", "",
        "Condition on eligible pair populations, annotation/lookup/scoring options and",
        "the retained full precomputed mean mu_P. Assume uniform without-replacement",
        "selection and independent strata within a method, with fixed per-pair numeric",
        "return status/value under normal successful execution. If G of M missing-lookup",
        "pairs return numeric values, R is hypergeometric for c requests. The expected",
        "native ratio is theta=w*mu_P+(1-w)*mu_G, w=E[k/(k+R)], not substitution of E[R].",
        "Count-tail inversion and conditional bounded-mean concentration are projected",
        "through this mixture. Alpha=0.05 is allocated across two components for each",
        "of eight methods; their joint rectangle supplies all 28 difference ranges",
        "without assuming cross-method independence. Score-dependent omissions are",
        "allowed, with no missing-at-random or biological pair/family independence.", "",
        "| Retained Method | Observed Native Z | Conditional 95% Expected-Theta Range |",
        "| --- | ---: | --- |"]
    for row in report["methods"]:
        lo,hi = row["interval"]["expected_native_mean_bounds"]
        lines.append("| %s | %.9f | [%.9f, %.9f] |" % (row["method"],row["observed_native_mean"],lo,hi))
    lines += ["", "Observed Z and expected theta are distinct: unknown G/mu_G preclude an exact",
        "theta point value. These ranges are not biological error bars centered on Z.",
        "Both OrthoHMM modes exceed full OrthoFinder on this conditional FAS target,",
        "but Proteinortho exceeds both; high sensitivity exceeds phylogenetic OrthoHMM.",
        "SonicParanoid versus OrthoMCL overlaps zero. All favorable and unfavorable",
        "contrasts are retained. This is not overall orthology superiority or a causal",
        "explanation, and earlier negative F1 comparisons remain unchanged.", "",
        "The prospective finite-population validation covers 15 cells, including four",
        "literal 9000-cap cases; its context-dependent-return control is rejected.",
        "Independent 80-digit recurrence readback checks native count endpoints, weights",
        "and all ranges/projections without producer/kernel/SciPy imports. Conditional",
        "design assumptions, unestablished historical parser/database bindings and",
        "missing historical RNG/pair identities remain limitations, not certificates.",
        "No biological generalization, other-endpoint uncertainty or readiness is admitted.",
        "[Frozen reporting plan](NATIVE_FAS_CONDITIONAL_REPORTING_PROTOCOL_20261010.md),",
        "[actual conditional result and limits](NATIVE_FAS_CONDITIONAL_RESULT_20261010.md),",
        "[all 28 contrasts](native_fas_conditional_report_20261010_v1/summary.md),",
        "[actual execution and independent readback](native_fas_conditional_execution_20261010_v1.json).", "", ""]
    qualification = """The conditional repeated-statistic analysis above adds a distinct sampling
target for the retained eight-method FAS comparison. It does not identify
full-eligible-population functional accuracy or new-family/clade uncertainty.
Earlier population-bound and stratum-weight diagnostics alone do not supply
confidence intervals; the conditional ranges additionally require the stated
uniform/fixed-return design assumptions and validation. Historical source
identity and biological/generalization gaps remain unresolved. The newer
native factorial's GO/EC/FAS uncertainty is not supplied by this older snapshot.

"""
    return "\n".join(lines), qualification


def manuscript(parent, report, execution):
    require(all(parent.count(anchor) == 1 for anchor in ANCHORS), "Ambiguous manuscript anchors")
    inserts = sections(report,execution)
    revised = parent
    for anchor,section in zip(ANCHORS,inserts):
        require(section not in parent, "Already integrated")
        revised = revised.replace(anchor,section+anchor,1)
    restored = revised
    for section in inserts:
        restored = restored.replace(section,"",1)
    require(restored == parent, "Unscoped parent alteration")
    return revised,inserts


def run(root, output, receipt):
    directory = Path(root).resolve()/"benchmark_tools/results"
    output,receipt = Path(output).absolute(),Path(receipt).absolute()
    require(output != receipt and output.suffix == ".md" and receipt.suffix == ".json" and
            output.parent.resolve() == receipt.parent.resolve() == directory, "Preserve relative links and output scope")
    for path in (output,receipt):
        if path.exists() or path.is_symlink(): raise FileExistsError(path)
    refs,docs = {},{}
    for key,(name,sha) in PINS.items():
        path = directory/name
        refs[key] = record(path)
        require(refs[key]["sha256"] == sha,"Changed integration input: " + key)
        text = path.read_text()
        docs[key] = text if key == "parent" else json.loads(text)
    revised,inserts = manuscript(docs["parent"],docs["report"],docs["execution"])
    require(all(record(ref["path"]) == ref for ref in refs.values()),"Input changed during integration")
    with output.open("x") as stream: stream.write(revised)
    result = dict(schema="conditional_fas_manuscript_integration_v1",inputs=refs,source=record(__file__),
                  manuscript=record(output),inserted_sections=list(inserts),parent_unchanged_except_insertions=True,
                  conditional_design_ranges_reported=True,new_inference_or_scoring=False,
                  unconditional_historical_interval_admission=False,biological_generalization_intervals=False,
                  publication_ready=False,manuscript_rendered=False,visual_reviewed=False)
    with receipt.open("x") as stream:
        json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False); stream.write("\n")
    return result
