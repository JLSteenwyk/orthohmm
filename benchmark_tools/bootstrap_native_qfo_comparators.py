"""Retrospective native-cell comparisons with fixed 48-endpoint SwissTrees scope."""

import argparse
import csv
import io
import json
from pathlib import Path

import numpy as np

from benchmark_tools.audit_qfo_swiss_counts import LABELS, statistics
from benchmark_tools.bootstrap_corrected_swiss_comparators import validated_values
from benchmark_tools.bootstrap_qfo_factorial import CELLS
from benchmark_tools.bootstrap_qfo_swiss_stages import METRICS, aggregate
from benchmark_tools.bundle_publication_package import identity, pin, save
from benchmark_tools.export_native_factorial_progress import load, require


COMPARATORS = ("orthofinder_3_1_5_full", "orthofinder_3_1_5_sequence_only")
REPLICATES, SEED, MULTIPLICITY = 100000, 20260920, 48
PROTOCOL = "benchmark_tools/results/NATIVE_QFO_COMPARATOR_UNCERTAINTY_PROTOCOL_20261009.md"
PROTOCOL_SHA = "8941ed3c9fa206a988271faccc62484f9b1e4d55ef8195488024a36fc5c6a6ff"
SOURCES = {
    "legacy": ("benchmark_tools/results/qfo_recovered_swiss_uncertainty_22178.json",
               "a599749d66433ec211ed1ab0e3a6a2a75eb99cb0dc57abc47c2d267930cbbe34"),
    "binding": ("benchmark_tools/results/native_qfo_four_cell_swiss_uncertainty_20261007_v1.json",
                "e8f554a178841da24a1892a7e42f8680f9fe8e4827dec6d6bfe4643cdb890bac"),
    "snapshot": ("benchmark_tools/results/native_qfo_terminal_failures_20261009_v1/report.json",
                 "863ebc334dbdb5a41d1a9e9bba56be82896f0ff26c98366047928c7ef802966f"),
}
HELPERS = {
    "audit_qfo_swiss_counts.py": "ff8e5314c98bdffce54732df94611e9d340d854e70392358c495b636b67452ae",
    "bootstrap_corrected_swiss_comparators.py": "a5c339960936d0a41a66023059e4390db46e0cdc505194d209f56af060b97050",
    "bootstrap_qfo_swiss_comparators.py": "c7b418c92a8b12dc06bb1a51880f7ef82c3451a8af4a876cd3c2cbd1802f9122",
    "bootstrap_qfo_swiss_stages.py": "0598773a8bd04846f24e4cb27ab7b3529410c11451cbb399b9148a320fd06435",
    "bootstrap_qfo_factorial.py": "be09876ae4de31b923818385bd3d40d8d5e7df177215fc45a2101c7d30d1c1fa",
    "export_native_factorial_progress.py": "af87cf899b8ae474822ea81a922643ab2c229275f99a1df57d5149b0d280448a",
    "bundle_publication_package.py": "948efb04abfa84aea1c6822d9145592a15d54416838245f0c1bf569f5d81bbab",
}
MISSING_REASONS = {
    "p0_c1_r1": "pre_native_cpu_binding_failure",
    "p1_c0_r0": "not_in_supplied_fresh_cell_snapshot; historical method remains separate",
    "p1_c1_r0": "native_success_but_scoring_OUT_OF_MEMORY_no_admission",
    "p1_c1_r1": "native_inference_SIGSEGV_no_admission",
}


def prepare(root):
    root = Path(root).resolve()
    evidence, values = [], {}
    for key, (path, sha) in SOURCES.items():
        values[key], _ = load(root / path, sha, evidence)
    protocol = root / PROTOCOL
    require(identity(protocol)["sha256"] == PROTOCOL_SHA, "Changed native comparator protocol")
    helpers = []
    for name, sha in HELPERS.items():
        path = Path(__file__).with_name(name)
        require(identity(path)["sha256"] == sha, "Changed numerical/provenance helper: " + name)
        helpers.append(dict(path=str(path), **identity(path)))
    legacy, binding, snapshot = (values[name] for name in ("legacy", "binding", "snapshot"))
    require(legacy["status"] == "corrected_swiss_comparison_intervals_audited"
            and legacy["complete_panel"] is True and legacy["multiplicity_endpoints"] == 24
            and legacy["protocol_controls_match"] is True, "Wrong retained comparator result")
    require(binding["schema"] == "allocated_native_qfo_retained_swiss_uncertainty_binding_v1"
            and binding["multiplicity_endpoints"] == 42
            and binding["new_bootstrap_draws"] == 0 and binding["publication_ready"] is False,
            "Wrong retained native count binding")
    require(snapshot["schema"] == "native_qfo_terminal_failure_reporting_v1"
            and snapshot["publication_ready"] is False
            and [row["index"] for row in snapshot["rows"]] == list(range(6, 13)), "Wrong current native snapshot")
    targets = {row["cell"]: row for row in snapshot["rows"] if row["accuracy_admitted"]}
    require(set(targets) == {"p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p1_c0_r1"}
            and set(binding["bound_cells"]) == set(targets), "Changed admitted cell inventory")
    for row in snapshot["rows"]:
        require(row["accuracy_admitted"] == all(v is not None for v in row["scores"].values())
                and len(row["scores"]) == 6, "Missing native score cannot be admitted")
    retained = legacy["reconstructed_counts"]
    native = {}
    for cell, bound in binding["bound_cells"].items():
        audit, ref = load(bound["count_audit"]["path"], bound["count_audit"]["sha256"], evidence)
        require(ref == bound["count_audit"] and audit["reference"] == retained["reference"]
                and audit["families"] == retained["families"]
                and audit["reference_relation_count"] == retained["reference_relation_count"]
                and all(audit[key] is False for key in (
                    "publication_ready", "independent_confirmation", "new_accuracy_or_resource_admission"))
                and audit["new_bootstrap_draws"] == 0, "Count audit scope/reference differs")
        matches = [row for row in audit["cells"] if row["cell"] == cell]
        require(len(matches) == 1, "Ambiguous native count record")
        row, target = matches[0], targets[cell]
        require(row["admission"] == target["admission"]
                and all(row[key] == target[key] == bound[key]
                        for key in ("index", "native_job_id"))
                and row["native_endpoint_f1"] == target["scores"]["SwissTrees"]
                and row["raw_file"] in audit["checked_inputs"], "Native count/admission binding differs")
        require(abs(row["aggregate"]["F1"] - row["native_endpoint_f1"]) <= 5e-8,
                "Count statistic differs from native endpoint")
        if cell == "p0_c0_r1":
            require(row["resources"] is None and row["timing_admitted"] is False
                    and row["timing_eligible"] is False, "Recovered timing was relabeled")
        native[cell] = row
    reference = retained["reference"]
    path = Path(reference["path"])
    require(path.resolve().is_relative_to(root) and identity(path) == pin(reference),
            "Changed SwissTrees reference bytes")
    evidence.append(reference)
    return native, retained, dict(evidence=evidence, helpers=helpers,
        protocol=dict(path=str(protocol), **identity(protocol)))


def native_values(native, retained):
    require(native and set(native).issubset(CELLS), "Require known available native cells")
    checked = validated_values(retained)
    require(all(name in checked for name in COMPARATORS), "Missing OrthoFinder comparator")
    families = retained["families"]
    reference = next(row for row in retained["methods"] if row["method"] == COMPARATORS[0])["families"]
    values = {name: checked[name] for name in COMPARATORS}
    for cell in CELLS:
        if cell not in native:
            continue
        row = native[cell]
        require(row["cell"] == cell and [r["family"] for r in row["families"]] == families,
                "Native family order/cell differs")
        current = []
        for fresh, original in zip(row["families"], reference):
            counts, old = fresh["counts_without_prior"], original["counts_without_prior"]
            require(set(counts) == set(LABELS) and all(type(v) is int and v >= 0 for v in counts.values()),
                    "Invalid native confusion counts")
            require(fresh["represented_genes"] == original["represented_genes"]
                    and counts["TP"] + counts["FN"] == old["TP"] + old["FN"]
                    and counts["FP"] + counts["TN"] == old["FP"] + old["TN"],
                    "Native reference members/truth differ")
            scores, stored = statistics(counts), fresh["statistics_with_prior"]
            require(set(stored) == set(METRICS) and all(type(stored[m]) in (int, float)
                    and np.isfinite(stored[m]) and abs(scores[m] - stored[m]) <= 1e-12 for m in METRICS),
                    "Native family statistic differs from counts")
            current.append([scores["PPV"], scores["TPR"]])
        array = np.asarray(current)
        point, stored = aggregate(array.mean(axis=0)), row["aggregate"]
        require(set(stored) == set(METRICS) and all(type(stored[m]) in (int, float)
                and np.isfinite(stored[m]) and abs(point[j] - stored[m]) <= 1e-12 for j, m in enumerate(METRICS)),
                "Native aggregate is not benchmark macro statistic")
        values[cell] = array
    return values


def analyze(native, retained, replicates=REPLICATES, seed=SEED):
    require(type(replicates) is int and replicates >= 100 and type(seed) is int and seed >= 0,
            "Invalid bootstrap controls")
    values = native_values(native, retained)
    n = len(retained["families"])
    multiplicities = np.random.Generator(np.random.PCG64(seed)).multinomial(
        n, np.full(n, 1 / n), size=replicates)
    draws = {name: aggregate(multiplicities @ value / n) for name, value in values.items()}
    points = {name: aggregate(value.mean(axis=0)) for name, value in values.items()}
    comparisons = []
    for cell in CELLS:
        for comparator in COMPARATORS:
            row = dict(candidate=cell, reference=comparator, metrics=None, family_differences=None,
                status="native_cell_unavailable", missing_reason=MISSING_REASONS.get(cell, "not_supplied"))
            if cell in native:
                delta = draws[cell] - draws[comparator]
                family_delta = aggregate(values[cell]) - aggregate(values[comparator])
                row.update(status="conditional_native_comparison_estimated", missing_reason=None,
                    metrics={metric: dict(difference=float(points[cell][j] - points[comparator][j]),
                        paired_percentile_ci=np.quantile(delta[:, j], [.025, .975], method="linear").tolist(),
                        bonferroni_percentile_ci=np.quantile(delta[:, j],
                            [.05 / (2 * MULTIPLICITY), 1 - .05 / (2 * MULTIPLICITY)], method="linear").tolist(),
                        family_wins=int(np.sum(family_delta[:, j] > 1e-10)),
                        family_ties=int(np.sum(np.abs(family_delta[:, j]) <= 1e-10)),
                        family_losses=int(np.sum(family_delta[:, j] < -1e-10)))
                        for j, metric in enumerate(METRICS)},
                    family_differences=[dict(family=family, **dict(zip(METRICS, delta.tolist())))
                        for family, delta in zip(retained["families"], family_delta)])
            comparisons.append(row)
    return dict(schema="native_qfo_comparator_sensitivity_v1", comparisons=comparisons,
        families=retained["families"], observed_cells=[cell for cell in CELLS if cell in native],
        family_values={name: value.tolist() for name, value in values.items()},
        point_estimates={name: dict(zip(METRICS, value.tolist())) for name, value in points.items()},
        planned_endpoints=MULTIPLICITY, estimated_endpoints=6 * len(native), replicates=replicates,
        seed=seed, alpha=.05, multiplicity_endpoints=MULTIPLICITY, quantile_method="linear",
        rng="numpy.PCG64 multinomial; shared family multiplicities", numpy_version=np.__version__,
        units="raw 0-to-1 differences", protocol_controls_match=replicates == REPLICATES and seed == SEED,
        new_bootstrap_draws=replicates, retained_intervals_reused=False, unobserved_cells_imputed=False,
        independent_confirmation=False, new_accuracy_or_resource_admission=False, publication_ready=False,
        limitations=[
            "Retrospective development-exposed conditional analysis, not selection-adjusted or independent confirmation.",
            "All 48 planned endpoints stay in this separate correction; earlier 24/42 analyses remain unchanged.",
            "Only 18 families; shared history/merged predictions can violate exchangeability and percentile coverage.",
            "Regenerated deterministic multiplicities are computation, not another independent dataset.",
            "Sequence-only OrthoFinder is an MCL-checkpoint diagnostic; R also changes prediction semantics.",
            "Intervals do not cover other QfO challenges, secondary mean, timings or missing native outcomes.",
            "No equivalence, general superiority, initial-HMM causal effect or new default follows."])


def tables(result):
    fields = ["candidate", "reference", "status", "metric", "difference", "nominal_low", "nominal_high",
              "adjusted_low", "adjusted_high", "family_wins", "family_ties", "family_losses"]
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    lines = ["# Native-Cell SwissTrees Comparator Sensitivity", "",
        "Retrospective conditional analysis: 18 families; 100000 shared draws, seed 20260920.",
        "Candidate minus reference; raw 0-to-1 units. All 48 planned endpoints adjusted; 24 estimated.",
        "Sequence-only OrthoFinder is an MCL-checkpoint diagnostic. No independent confirmation or default change.", "",
        "| Native Cell | Comparator | Metric | Difference | Nominal 95% CI | Adjusted CI | Wins/Ties/Losses |",
        "| --- | --- | --- | ---: | --- | --- | --- |"]
    for row in result["comparisons"]:
        for metric in METRICS:
            out = dict.fromkeys(fields, "NA")
            out.update({key: row[key] for key in ("candidate", "reference", "status")}, metric=metric)
            if row["metrics"] is not None:
                value = row["metrics"][metric]
                out.update(difference=value["difference"], nominal_low=value["paired_percentile_ci"][0],
                    nominal_high=value["paired_percentile_ci"][1], adjusted_low=value["bonferroni_percentile_ci"][0],
                    adjusted_high=value["bonferroni_percentile_ci"][1],
                    **{key: value[key] for key in ("family_wins", "family_ties", "family_losses")})
                cells = [f"{value['difference']:+.6f}",
                         *("[%.6f, %.6f]" % tuple(value[key]) for key in (
                             "paired_percentile_ci", "bonferroni_percentile_ci")),
                         f"{value['family_wins']}/{value['family_ties']}/{value['family_losses']}"]
            else:
                cells = ["NA", "NA", "NA", "NA"]
            writer.writerow(out)
            lines.append("| " + " | ".join([row["candidate"], row["reference"], metric, *cells]) + " |")
    lines.extend(["", *["- " + note for note in result["limitations"]], ""])
    return stream.getvalue(), "\n".join(lines)


def run(root, output):
    output = Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    native, retained, provenance = prepare(root)
    result = analyze(native, retained)
    require(result["protocol_controls_match"] and result["estimated_endpoints"] == 24,
            "Production protocol controls differ")
    result.update(provenance, source=dict(path=str(Path(__file__).resolve()), **identity(Path(__file__))))
    for ref in [*provenance["evidence"], *provenance["helpers"], provenance["protocol"], result["source"]]:
        require(identity(ref["path"]) == pin(ref), "Input/helper changed during analysis")
    tsv, markdown = tables(result)
    output.mkdir(parents=True)
    save(output / "report.json", result)
    (output / "intervals.tsv").write_text(tsv)
    (output / "TABLE.md").write_text(markdown)
    return dict(report=identity(output / "report.json"), planned_endpoints=MULTIPLICITY,
                estimated_endpoints=24, publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.root, args.output)))
