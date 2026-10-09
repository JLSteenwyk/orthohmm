"""Insert verified controlled-fragment evidence without changing the frozen draft."""

import argparse
import json
import math
from pathlib import Path

from benchmark_tools.export_native_qfo_three_cell_strata import record, require, write_tsv


PINS = {
    "result": ("benchmark_tools/results/controlled_fragment_results_20261009_v1/report.json",
               "b1480e8eb2ace2aed743f5fc4e515e27ac4ba68b48e644b33989eda3cf8a7ab7"),
    "reader": ("benchmark_tools/results/controlled_fragment_execution_20261009_v1.json",
               "af665df1fff1ed454c9cce0286d9bcac831e935ad50b65aba5ee14e76d58edcf"),
    "panel": ("benchmarks/work/controlled_fragment_observations_20261009_v1/manifest.json",
              "8e58830911f60de5aec9ad616e5c02c12229702780e26216aa1371017637666d"),
    "runtime": ("benchmark_tools/results/controlled_fragment_execution_runtime_20261009_v1.json",
                "4ea5d10bccd02ef42e49194908212d07af701c109db878e176f56c298c2ca4da"),
    "parent": ("benchmark_tools/results/native_qfo_four_cell_strata_manuscript_20261009_v1.md",
               "60af8bf7ead43a3860eee6c923fc793895a8d8ae8e5564a30b0fe0daef131efc"),
}
SEEDS = list(range(20261101, 20261111))
METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2",
           "orthofinder_full", "orthofinder_sequence_only")
LABELS = dict(zip(METHODS, ("HS", "PHY", "OF", "OF checkpoint")))
METRICS = ("f1", "precision", "recall")
ARMS = ("baseline", "fragment")
ANCHORS = ("### Fresh Native Ablation Protocol\n",
           "### Gene-Tree And Event-History Controls Localize Errors\n",
           "## Reproducibility And Availability\n")
MEAN_FIELDS = ("method", "baseline_f1", "fragment_f1", "fragment_precision", "fragment_recall",
               "fragment_pair_endpoint_coverage", "eligible_seeds")
DIFF_FIELDS = ("target_arm", "target_method", "reference_arm", "reference_method", "metric",
               "estimate", "nominal_low", "nominal_high", "adjusted_low", "adjusted_high", "eligible_seeds")
COUNT_FIELDS = ("method", "fragment_endpoints", "baseline_tp", "baseline_fp", "baseline_fn",
                "fragment_tp", "fragment_fp", "fragment_fn")


def bind(ref, actual):
    return (ref["absolute_path"] == actual["path"] and ref["bytes"] == actual["bytes"]
            and ref["sha256"] == actual["sha256"])


def validated_tables(docs):
    result, reader, panel, runtime = (docs[k] for k in ("result", "reader", "panel", "runtime"))
    require(result["schema"] == "controlled_fragment_results_v1"
            and result["status"] == "complete_with_explicit_outcomes"
            and result["publication_ready"] is False and result["fragment_successes"] == 40,
            "Incomplete fragment evidence")
    require(reader["schema"] == "controlled_fragment_actual_execution_v1"
            and reader["publication_ready"] is False and reader["automatic_inference_retries"] == 0,
            "Changed executed scope")
    terminal = reader["independent_result_readback"]["terminal"]
    require(terminal["exit_code"] == 0, "Independent reader did not succeed")
    checked = json.loads(terminal["output"])
    require(checked["status"] == "independent_counts_intervals_and_tables_verified"
            and [checked[k] for k in ("score_records", "stratum_records", "planned_comparisons",
                                      "bootstrap_replicates")] == [80, 240, 15, 20000]
            and checked["report_sha256"] == PINS["result"][1], "Incomplete independent readback")
    require(panel["schema"] == "controlled_fragment_observations_v1"
            and panel["condition"] == "fragment20_center60_v1" and panel["inference_identities"] == 30
            and [d["seed"] for d in panel["datasets"]] == SEEDS, "Changed fragment panel")
    amendment = runtime["inventory_amendment"]
    require(runtime["schema"] == "controlled_fragment_execution_runtime_v1"
            and amendment["historical_output_equivalence_established"] is False
            and amendment["historical_full_inventory_equal"] is False
            and len(amendment["differences"]) == 72
            and len(amendment["scientific_dependencies_unchanged"]) == 11, "Changed runtime scope")
    job = result["scheduler"]
    require((job["JobIDRaw"], job["State"], job["ExitCode"], job["AllocCPUS"], job["ReqMem"])
            == ("24087", "COMPLETED", "0:0", "16", "16G"), "Changed terminal execution")
    records = {(r["arm"], r["seed"], r["method"]): r for r in result["records"]}
    require(len(records) == len(result["records"]) == 80
            and set(records) == {(a, s, m) for a in ARMS for s in SEEDS for m in METHODS},
            "Missing or duplicate result identities")
    for row in records.values():
        require(row["status"] == "complete" and [s["fragment_endpoints"] for s in row["strata"]] == [0, 1, 2],
                "Unavailable or changed endpoint strata")
        for score in (row["score"], *(s["score"] for s in row["strata"])):
            require(all(type(score[k]) is int and score[k] >= 0 for k in ("tp", "fp", "fn"))
                    and all(type(score[k]) in (int, float) and math.isfinite(score[k])
                            and 0 <= score[k] <= 1 for k in (*METRICS, "pair_endpoint_coverage")),
                    "Invalid complete score")
        require(all(sum(s["score"][k] for s in row["strata"]) == row["score"][k]
                    for k in ("tp", "fp", "fn")), "Endpoint strata do not partition counts")
    summaries = {(r["arm"], r["method"]): r for r in result["metric_means"]}
    require(len(summaries) == len(result["metric_means"]) == 8
            and set(summaries) == {(a, m) for a in ARMS for m in METHODS}, "Incomplete method means")
    for (arm, method), row in summaries.items():
        require(row["planned_seeds"] == 10 and row["successful_seeds"] == SEEDS
                and row["failed_or_unavailable_seeds"] == [], "Changed summary eligibility")
        for metric in METRICS:
            expected = sum(records[arm, s, method]["score"][metric] for s in SEEDS) / 10
            require(row["metrics"][metric]["eligible_seeds"] == SEEDS
                    and math.isclose(row["metrics"][metric]["mean"], expected, rel_tol=0, abs_tol=1e-15),
                    "Changed mean or seed list")
    planned = [("fragment", m, "baseline", m) for m in METHODS[:3]] + [
        ("fragment", m, "fragment", METHODS[2]) for m in METHODS[:2]]
    comparisons = result["comparisons"]
    require([(r["target_arm"], r["target_method"], r["reference_arm"], r["reference_method"], r["metric"])
             for r in comparisons] == [(*p, metric) for p in planned for metric in METRICS],
            "Missing or changed planned comparisons")
    differences = []
    for row in comparisons:
        require(row["eligible_seeds"] == SEEDS and row["excluded_seeds"] == []
                and row["status"] == "conditional_approximation" and row["planned_endpoints"] == 15
                and row["replicates"] == 20000 and row["rng"] == "PCG64"
                and row["rng_seed"] == 20261011 and row["quantile_method"] == "linear",
                "Changed interval scope")
        values = [records[row["target_arm"], s, row["target_method"]]["score"][row["metric"]]
                  - records[row["reference_arm"], s, row["reference_method"]]["score"][row["metric"]]
                  for s in SEEDS]
        require(row["paired_differences"] == values
                and math.isclose(row["estimate"], sum(values) / 10, rel_tol=0, abs_tol=1e-15),
                "Changed paired effect")
        value = {k: row[k] for k in DIFF_FIELDS[:6]}
        for kind in ("nominal", "adjusted"):
            interval = row[kind + "_interval"]
            require(len(interval) == 2 and all(type(v) in (int, float) and math.isfinite(v)
                    and -1 <= v <= 1 for v in interval) and interval[0] <= interval[1], "Invalid interval")
            value.update({kind + "_low": interval[0], kind + "_high": interval[1]})
        require(value["adjusted_low"] <= value["nominal_low"] <= value["nominal_high"]
                <= value["adjusted_high"], "Non-nested intervals")
        differences.append(dict(value, eligible_seeds=10))
    for row in differences:
        if row["metric"] == "f1":
            require((row["adjusted_low"] <= 0 <= row["adjusted_high"] if row["reference_arm"] == "baseline"
                     else row["estimate"] < 0 and row["adjusted_high"] < 0), "Changed F1 interpretation")
    means, counts = [], []
    for method in METHODS:
        means.append(dict(method=method, baseline_f1=summaries["baseline", method]["metrics"]["f1"]["mean"],
            **{"fragment_" + m: summaries["fragment", method]["metrics"][m]["mean"] for m in METRICS},
            fragment_pair_endpoint_coverage=sum(records["fragment", s, method]["score"]["pair_endpoint_coverage"]
                                                for s in SEEDS) / 10, eligible_seeds=10))
        for endpoints in range(3):
            counts.append(dict(method=method, fragment_endpoints=endpoints,
                **{a + "_" + k: sum(records[a, s, method]["strata"][endpoints]["score"][k] for s in SEEDS)
                   for a in ARMS for k in ("tp", "fp", "fn")}))
    return means, differences, counts


def sections(docs):
    means, differences, counts = validated_tables(docs)
    methods = """### Controlled Fragment Observation Test

One prospectively frozen observation condition reuses the ten variable-length
baseline histories (seeds 20261101--20261110), not new biological histories.
Within each dataset, SHA256 of UTF-8 `seed:condition:gene_id` ranks all gene IDs;
the first floor(N/5), with ID tie-breaking, retain a centered floor(3L/5)
substring starting at floor((L-floor(3L/5))/2). Identifiers, species ownership,
unselected sequences and the byte-exact parent orthology truth stay unchanged.
Independent readback verified all 8,200 genes and 1,636 truncations. This is
controlled synthetic truncation, not natural fragment annotation or a new
independent generalization set. The simulator does not add within-family indels
or validated domain architectures in this condition.

All 30 fresh native runs completed: high-sensitivity OrthoHMM (HS), satellite_v2
phylogenetic OrthoHMM (PHY), and full OrthoFinder 3.1.5 (OF), sequentially for
each seed. The OF sequence-only MCL checkpoint is an additional diagnostic
output, not a separately rerun final sequence-only workflow or timing estimate.
Forty retained baseline outcomes are reused without rerunning controls. Settings
are frozen except input/output/metrics paths. The OrthoHMM CPU budget is four
with four threads per worker; actual metadata reports one search worker and
four total search threads, not four simultaneous workers. The sequential job
received 16 CPUs and 16 GiB on the shared Threadripper. The prospective runtime
amendment records 72 full-inventory differences, while eleven required
scientific dependency versions and the full OF inventory remain unchanged.
This is not proof of historical OrthoHMM output equivalence.

Micro cross-species pair F1, precision, recall and pair-endpoint coverage are
calculated per seed; displayed means are arithmetic means across the ten
seeds, not pooled-pair F1 or a harmonic mean of mean precision and recall.
Coverage is the fraction of all input genes present in at least one predicted
cross-species pair, not recall. Zero/one/two truncated-endpoint strata use the
same prospective flags in both arms and include inter-origin-family false
positives. Their integer TP/FP/FN counts partition each whole-dataset count.

The fixed 15-endpoint comparison family comprises fragment-minus-baseline for
HS, PHY and OF plus fragment HS-minus-OF and PHY-minus-OF, each for F1,
precision and recall. Whole seeds, not dependent gene pairs, are resampled
20,000 times with PCG64 seed 20261011 and linear quantiles. Nominal intervals
use .025/.975 and fixed-15 Bonferroni intervals use .05/30 and 1-.05/30.
These are conditional percentile approximations with only ten seed units;
extreme-tail coverage is not guaranteed. All planned seed lists are complete.
No checkpoint comparison interval or subgroup interval is added after seeing
results. [Frozen protocol](CONTROLLED_FRAGMENT_OBSERVATION_PROTOCOL_20261009.md),
[runtime amendment](CONTROLLED_FRAGMENT_RUNTIME_AMENDMENT_20261009.md),
[executed validation](controlled_fragment_execution_20261009_v1.json).

"""
    lines = ["### Controlled Fragment Accuracy", "",
        "All 40 fragment method/checkpoint outcomes were admitted across ten seeds.",
        "HS denotes high sensitivity, PHY denotes satellite_v2 with phylogeny, and",
        "OF denotes full OrthoFinder. Values below are percent; coverage is defined",
        "in the controlled observation protocol above.", "",
        "| Method | Baseline F1 | Fragment F1 | Fragment Precision | Fragment Recall | Coverage |",
        "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for row in means:
        lines.append("| " + " | ".join([LABELS[row["method"]], *[f"{100*row[k]:.4f}" for k in MEAN_FIELDS[1:-1]]]) + " |")
    lines.extend(["", "### Controlled Fragment Paired Differences", "",
        "All 15 planned differences and interval endpoints are shown in percentage",
        "points. F-B denotes fragment minus baseline; F-F denotes the fragment",
        "method minus fragment OF. Each row has ten eligible paired seeds.", "",
        "| Comparison | Metric | Difference | Nominal 95% Interval | Fixed-15 Adjusted Interval | Seeds |",
        "| --- | --- | ---: | --- | --- | ---: |"])
    for row in differences:
        comparison = (LABELS[row["target_method"]] + " F-B" if row["reference_arm"] == "baseline"
                      else LABELS[row["target_method"]] + "-OF F-F")
        lines.append("| " + " | ".join([comparison, row["metric"], f"{100*row['estimate']:+.4f}",
            *[f"[{100*row[k+'_low']:+.4f}, {100*row[k+'_high']:+.4f}]" for k in ("nominal", "adjusted")], "10"]) + " |")
    lines.extend(["", "Both fragment HMM configurations have lower mean F1 than full OF in this",
        "condition; both adjusted F1 difference intervals exclude zero on the",
        "negative side. The three native fragment-minus-baseline F1 intervals",
        "include zero, which establishes neither equivalence nor absence of",
        "fragment sensitivity. No default was changed based on these outcomes.", "",
        "### Controlled Fragment Endpoint Counts", "",
        "Counts below sum the ten seeds; they are not the mean-score or interval",
        "statistics. Endpoint count is the number of prospectively flagged genes",
        "in a pair, including the reused untruncated baseline arm.", "",
        "| Method | Truncated Endpoints | Baseline TP / FP / FN | Fragment TP / FP / FN |",
        "| --- | ---: | --- | --- |"])
    for row in counts:
        lines.append("| " + " | ".join([LABELS[row["method"]], str(row["fragment_endpoints"]),
            *[" / ".join(str(row[a + "_" + k]) for k in ("tp", "fp", "fn")) for a in ARMS]]) + " |")
    lines.extend(["", "Unflagged-pair counts can change through grouping and reconciliation even",
        "though neither endpoint was truncated. These counts alone do not localize",
        "a causal stage. Representative search/candidate/group/reconciliation tracing",
        "remains separate work. The initial HMM search is on in both HMM methods;",
        "this observation test is not an initial-HMM-off causal comparison.",
        "[All 80 scores, 240 strata and 15 comparisons](controlled_fragment_results_20261009_v1/report.json),",
        "[full-precision stratum table](controlled_fragment_results_20261009_v1/strata.tsv),",
        "[result interpretation](CONTROLLED_FRAGMENT_RESULT_20261009.md).", "", ""])
    limitations = """### Controlled Fragment Scope And Limitations

The controlled observation result supports one bounded synthetic-fragment
diagnostic, not a claim of better natural-fragment recovery or broad superiority
to full OrthoFinder. It preserves the negative comparisons. Development
exposure, artificial centered truncation, absence of validated domain/indel
truth, ten-seed approximate interval tails and the changed full OrthoHMM
runtime inventory limit interpretation. Counts that change in the unflagged
stratum do not identify which search, grouping or reconciliation step caused
the error. No initial-HMM-off control is supplied by this extension.

The fragment job is not a replacement scaling panel. Cached baseline costs
are not paired timings, and the checkpoint has no independent runtime.
Timing measurements were collected on a shared Threadripper while other
analyses were running. Competition for CPU, memory bandwidth and I/O may have
affected elapsed times, with an unknown and potentially tool-dependent impact.
These are observed shared-host timings, not estimates of isolated performance.
This extension does not resolve unavailable original TreeFam resources,
failed native factorial cells or other unmet publication requirements.

"""
    return (methods, "\n".join(lines), limitations), (means, differences, counts)


def manuscript(parent, docs):
    require(all(parent.count(anchor) == 1 for anchor in ANCHORS), "Ambiguous manuscript anchors")
    inserts, tables = sections(docs)
    require(all(section not in parent for section in inserts), "Already integrated")
    revised = parent
    for anchor, section in zip(ANCHORS, inserts):
        revised = revised.replace(anchor, section + anchor, 1)
    restored = revised
    for section in inserts:
        restored = restored.replace(section, "", 1)
    require(restored == parent, "Changed frozen parent body")
    return revised, inserts, tables


def run(root, output, revised_path):
    root = Path(root).resolve()
    directory = root / "benchmark_tools/results"
    output, revised_path = Path(output).absolute(), Path(revised_path).absolute()
    for path in (output, revised_path):
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    require(output.parent.resolve() == revised_path.parent.resolve() == directory, "Preserve manuscript-relative links")
    docs, refs = {}, {}
    for key, (name, sha) in PINS.items():
        path = root / name
        refs[key] = record(path)
        require(refs[key]["sha256"] == sha, "Changed integrated input: " + key)
        docs[key] = path.read_text() if key == "parent" else json.loads(path.read_text())
    require(bind(docs["result"]["manifest"], refs["panel"])
            and bind(docs["result"]["runtime"], refs["runtime"])
            and docs["reader"]["outputs"]["sha256"] == refs["result"]["sha256"]
            and docs["reader"]["prepared_manifest"]["sha256"] == refs["panel"]["sha256"]
            and docs["reader"]["runtime_manifest"]["sha256"] == refs["runtime"]["sha256"],
            "Mixed executed evidence bindings")
    revised, inserts, tables = manuscript(docs["parent"], docs)
    source = record(__file__)
    for ref in [*refs.values(), source]:
        require(record(ref["path"]) == ref, "Changed input before writing")
    output.mkdir(parents=True)
    for name, section in zip(("methods", "results", "limitations"), inserts):
        (output / (name + "_section.md")).write_text(section)
    for name, rows, fields in zip(("means", "comparisons", "endpoint_counts"), tables,
                                 (MEAN_FIELDS, DIFF_FIELDS, COUNT_FIELDS)):
        write_tsv(output / (name + ".tsv"), rows, fields)
    (output / "claims.md").write_text("""# Controlled Fragment Claim Addendum

- Supported: one frozen synthetic observation condition completed all 30 native
  processes; 40 fragment method/checkpoint outcomes and all controls were verified.
- Supported: full OF has higher mean F1 than both HMM configurations in this
  condition, including negative adjusted difference intervals. No superiority claim.
- Unestablished: zero-containing fragment-minus-baseline F1 intervals are not
  equivalence or proof of insensitivity; no subgroup or checkpoint intervals.
- Unestablished: natural fragment recovery, complete domain architecture truth,
  initial-HMM causal benefit, isolated efficiency or general publication readiness.
- Evidence: ../controlled_fragment_results_20261009_v1/report.json,
  ../controlled_fragment_execution_20261009_v1.json and this generated Methods,
  Results, limitations, mean/comparison/endpoint-count tables. Earlier dated
  claim checklists remain unchanged; representative pipeline tracing is unfinished.
""")
    revised_path.write_text(revised)
    for ref in [*refs.values(), source]:
        require(record(ref["path"]) == ref, "Changed input while writing")
    manifest = dict(schema="controlled_fragment_integration_v1", inputs=refs, source=source,
        outputs=[record(p) for p in sorted(output.iterdir())], manuscript=record(revised_path),
        table_rows=dict(means=4, comparisons=15, endpoint_counts=12),
        parent_body_unchanged_except_insertions=True, table_scale=100, tsv_scale=1,
        new_bootstrap_draws=0, new_inference_or_scoring=False, publication_ready=False,
        manuscript_rendered=False, visual_reviewed=False)
    with (output / "manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output", "manuscript"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.root, args.output, args.manuscript), allow_nan=False))
