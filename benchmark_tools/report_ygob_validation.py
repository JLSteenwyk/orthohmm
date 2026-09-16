"""Report assembly for the frozen YGOB protocol; callers must verify run gates."""

from benchmark_tools.orthofinder_mcl_to_orthogroups import iter_mcl_clusters, load_sequence_ids
from benchmark_tools.score_ygob_groups import membership, paired_bootstrap, score_groups


METHODS = {
    "orthohmm_satellite_v2": "OrthoHMM satellite_v2",
    "orthohmm_high_sensitivity": "OrthoHMM high sensitivity",
    "orthofinder_full": "OrthoFinder 3.1.5 full",
    "orthofinder_sequence_only": "OrthoFinder sequence-only checkpoint (diagnostic)",
}
DIAGNOSTIC = "orthofinder_sequence_only"


def read_checkpoint(clusters_path, sequence_ids_path, input_genes):
    """Require one-to-one ID restoration and complete checkpoint membership."""
    mapping = load_sequence_ids(sequence_ids_path)
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("Multiple internal IDs map to the same original gene")
    if set(mapping.values()) != set(input_genes):
        raise ValueError("SequenceIDs must match the complete inference universe")
    groups = {}
    for index, cluster in enumerate(iter_mcl_clusters(clusters_path)):
        if any(gene not in mapping for gene in cluster):
            raise ValueError("Unknown internal ID in MCL checkpoint")
        groups[f"MCL{index}"] = [mapping[gene] for gene in cluster]
    if set(membership(groups)) != set(input_genes):
        raise ValueError("MCL checkpoint does not cover the complete inference universe")
    return groups


def assemble_report(predictions, references, input_genes):
    """Compute fixed endpoints, without selecting methods or tuning parameters."""
    if set(predictions) != set(METHODS):
        raise ValueError("Exactly the four frozen methods are required")
    input_genes = set(input_genes)
    scores = {method: score_groups(predictions[method], references, input_genes)
              for method in METHODS}
    uncertainty = paired_bootstrap({m: s for m, s in scores.items() if m != DIAGNOSTIC},
                                   baseline="orthofinder_full", replicates=20000,
                                   seed=20260917, batch_size=128)
    if uncertainty["multiplicity_count"] != 6:
        raise AssertionError("Frozen multiplicity changed")
    return {
        "schema_version": 1,
        "protocol": "benchmark_tools/results/YGOB_VALIDATION_PROTOCOL_20260916.md",
        "publication_ready": False,
        "completion_gates_verified_by_this_module": False,
        "input_gene_count": len(input_genes),
        "scores": scores, "uncertainty": uncertainty,
        "primary_contrast": "orthohmm_satellite_v2 versus orthofinder_full",
        "secondary_contrast": "orthohmm_high_sensitivity versus orthofinder_full",
        "diagnostic_methods": [DIAGNOSTIC],
        "limitations": [
            "Curated homolog-group co-membership, including within-species pairs, not resolved pairwise orthology.",
            "Novel-taxon transfer does not establish family-disjoint or unrestricted generalization.",
            "Predictions outside the reference universe are projected out, not counted as false positives.",
            "Coverage and missing predictions must be interpreted alongside F1.",
            "Run completion, frozen inputs and versions, overlap and resource audits require separate verification.",
            *uncertainty["limitations"],
        ],
    }


def markdown(report):
    lines = ["# YGOB Curated Group Recovery", "",
             "Completion gates require separate verification; this report alone does not establish publication readiness.", "",
             "| Method | F1 (%) | Precision (%) | Recall (%) | Reference-gene coverage (%) | Predicted-group coverage (%) | Exact reference groups |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for method, label in METHODS.items():
        score = report["scores"][method]
        values = [100 * score["metrics"][m] for m in ("f1", "precision", "recall")]
        values += [100 * score["reference_gene_coverage"], 100 * score["predicted_group_coverage"]]
        lines.append(f"| {label} | " + " | ".join(f"{v:.6f}" for v in values)
                     + f" | {score['exact_reference_groups']}/{score['reference_groups']} |")
    lines += ["", "Predicted-group coverage is the fraction of supplied groups retaining at least one scored gene.",
              "", "## Paired Differences Versus Full OrthoFinder", "",
              "20,000 paired pillar replicates; PCG64 seed 20260917. Differences and intervals are percentage points.", "",
              "| Contrast | Metric | Difference | Nominal 95% CI | Bonferroni CI (six contrasts/metrics) |",
              "| --- | --- | ---: | --- | --- |"]
    for method in ("orthohmm_satellite_v2", "orthohmm_high_sensitivity"):
        role = "Primary" if method == "orthohmm_satellite_v2" else "Secondary"
        for metric, result in report["uncertainty"]["comparisons"][method].items():
            intervals = [", ".join(f"{v:.6f}" for v in result[key])
                         for key in ("paired_95_percent_ci", "bonferroni_ci")]
            lines.append(f"| {role}: {METHODS[method]} | {metric} | "
                         f"{result['difference_percentage_points']:.6f} | [{intervals[0]}] | [{intervals[1]}] |")
    lines += ["", "## Limitations", "", *["- " + note for note in report["limitations"]]]
    return "\n".join(lines) + "\n"
