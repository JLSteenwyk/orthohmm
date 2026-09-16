#!/usr/bin/env python3
"""Consolidate retained benchmark evidence without treating missing runs as scores."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.qfo_summarize_scores import METRICS, metric_path, summarize


# Stable identifiers and explicit crosswalks, never inferred from rankings.
METHODS = (
    ("orthohmm_high_sensitivity", "OrthoHMM high sensitivity", "orthohmm_high_sensitivity", "group-derived pairs"),
    ("orthohmm_phylogeny_satellite_v2", "OrthoHMM phylogeny satellite_v2", "orthohmm_phylogeny_satellite_v2", "phylogenetically inferred pairs"),
    ("orthofinder_3_1_5_full", "OrthoFinder 3.1.5 full", "orthofinder_v3_diamond", "phylogenetically inferred pairs"),
    ("orthofinder_3_1_5_sequence_only", "OrthoFinder 3.1.5 sequence-only checkpoint", "orthofinder_v3_sequence_only", "MCL checkpoint group-derived pairs"),
    ("sonicparanoid_2_0_9", "SonicParanoid 2.0.9", "sonicparanoid_native", "native species-pair relations"),
    ("proteinortho_6_3_6", "ProteinOrtho 6.3.6", "proteinortho_native", "native post-clustering graph"),
    ("fastoma_0_3_5", "FastOMA 0.3.5 final orthologous groups", "fastoma", "native pairs; supplied OrthoFinder species tree"),
    ("orthomcl_1_4", "OrthoMCL 1.4", "orthomcl_1_4_final_groups", "final MCL group-derived pairs"),
)


def unique_index(records, key):
    result = {}
    for record in records:
        name = record[key]
        if name in result:
            raise ValueError(f"Duplicate {key}: {name}")
        result[name] = record
    return result


def qfo_record(directory):
    missing = [relative for relative, _ in METRICS.values()
               if not (directory / "results" / relative).is_file()]
    if missing:
        return {"status": "incomplete", "missing_metrics": missing}
    summary = summarize(directory)
    if any(not math.isfinite(value) or not 0 <= value <= 1 for value in summary["scores"].values()):
        raise ValueError(f"Out-of-range QfO score in {directory}")
    details = {}
    for name, (relative, _) in METRICS.items():
        path = metric_path(directory, relative)
        data = json.loads(path.read_text())["datalink"]["inline_data"]
        participants = data["challenge_participants"]
        if len(participants) != 1 or participants[0]["participant_id"] != directory.name:
            raise ValueError(f"Unexpected participant in {path}")
        for axis in ("metric_x", "metric_y"):
            value = participants[0][axis]
            if not math.isfinite(value) or value < 0 or (relative.startswith(("VGNC/", "SwissTrees/", "TreeFam-A/")) and value > 1):
                raise ValueError(f"Invalid QfO axis {axis} in {path}")
        details[name] = {"participant": participants[0],
                         "axes": data["visualization"], "source": file_provenance(path)}
    return {"status": "metrics_available", **summary, "metric_details": details}


def validate_three_score(score):
    tp, fp, fn = [score[k] for k in ("true_positive_gene_pairs", "false_positive_gene_pairs", "false_negative_gene_pairs")]
    if any(not isinstance(v, int) or v < 0 for v in (tp, fp, fn)):
        raise ValueError("Invalid Three Kingdoms counts")
    expected = {"precision": tp / (tp + fp) if tp + fp else 0,
                "recall": tp / (tp + fn) if tp + fn else 0,
                "f_score": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0}
    for key, value in expected.items():
        if not math.isclose(score[key], value, abs_tol=1e-12, rel_tol=0):
            raise ValueError(f"Three Kingdoms {key} disagrees with counts")


def verify_final_group_workflow(directory):
    """Require successful conversion/scoring and verify both submitted pair files."""
    manifest = directory / "result.sha256"
    if not manifest.is_file():
        return None
    metadata_path = directory / "run_metadata.tsv"
    metadata = {}
    for line in metadata_path.read_text().splitlines():
        key, value = line.split("\t", 1)
        if key in metadata:
            raise ValueError(f"Duplicate workflow metadata: {key}")
        metadata[key] = value
    if metadata.get("exit_code") != "0" or not metadata.get("finished"):
        raise ValueError("Final-group workflow did not finish successfully")
    expected = {(directory / name).resolve() for name in ("pairs.tsv", "pairs.qfo.tsv")}
    verified = {}
    for line in manifest.read_text().splitlines():
        digest, filename = line.split("  ", 1)
        path = Path(filename).resolve()
        if path not in expected or str(path) in verified:
            raise ValueError(f"Unexpected or duplicate final-group artifact: {path}")
        provenance = file_provenance(path)
        if provenance["sha256"] != digest:
            raise ValueError(f"Final-group artifact checksum mismatch: {path}")
        verified[str(path)] = provenance
    if set(map(Path, verified)) != expected:
        raise ValueError("Incomplete final-group artifact manifest")
    return {"metadata": metadata, "metadata_source": file_provenance(metadata_path),
            "manifest": file_provenance(manifest), "artifacts": list(verified.values())}


def build_report(root):
    results = root / "benchmark_tools/results"
    comparison_path = results / "external_tool_comparison_20260904.json"
    three_path = results / "three_kingdoms_parity_20260907.json"
    uncertainty_path = results / "orthobench_paired_uncertainty_20260916.json"
    historical = json.loads(comparison_path.read_text())
    three = json.loads(three_path.read_text())
    uncertainty = json.loads(uncertainty_path.read_text())
    ob_index = unique_index(historical["orthobench"]["records"], "method")
    three_index = unique_index(three["methods"], "key")
    records = []
    for key, label, qfo_dir, semantics in METHODS:
        ob = dict(ob_index[label])
        if key in uncertainty["scores"]:
            current = uncertainty["scores"][key]
            for metric in ("f_score", "precision", "recall"):
                if not math.isclose(ob[metric + "_percent"], current[metric], abs_tol=1e-8):
                    raise ValueError(f"OrthoBench baseline changed: {key} {metric}")
                ob[metric + "_percent"] = current[metric]
            ob["prediction_provenance"] = uncertainty["inputs"]["predictions"][key]
            ob["uncertainty"] = uncertainty["comparisons"].get(key)
        validate_three_score(three_index[key]["score"])
        qfo = qfo_record(root / "qfo_benchmark/scoring" / qfo_dir)
        qfo["output_semantics"] = semantics
        if key == "orthomcl_1_4":
            completion = verify_final_group_workflow(root / "qfo_benchmark/results" / qfo_dir)
            if completion is None:
                qfo["status"] = "awaiting_workflow_completion"
            else:
                qfo["workflow_completion"] = completion
        records.append({"key": key, "method": label, "qfo": qfo,
                        "orthobench": ob, "three_kingdoms": three_index[key]})
    diagnostic = qfo_record(root / "qfo_benchmark/scoring/orthomcl_1_4")
    diagnostic.update(output_semantics="pre-MCL cross-species weight-matrix edges", primary=False)
    return {
        "schema_version": 1, "generated_at": datetime.now(timezone.utc).isoformat(),
        "publication_ready": False, "source": file_provenance(Path(__file__)),
        "command": [sys.executable, *sys.argv],
        "inputs": [file_provenance(p) for p in (comparison_path, three_path, uncertainty_path)],
        "methods": records, "diagnostics": {"orthomcl_preclustering": diagnostic},
        "output_semantics_sources": {
            "proteinortho_6_3_6": "https://gitlab.com/paulklemm_PHD/proteinortho/-/raw/v6.3.6/README.md",
            "orthofinder_3_1_5_full": "benchmark_tools/results/ORTHOBENCH_UNCERTAINTY_PROTOCOL_20260916.md",
            "orthomcl_1_4": "benchmark_tools/results/PUBLICATION_PROGRESS.md",
        },
        "limitations": [
            "QfO and OrthoBench are development-exposed; these results do not establish independent generalization.",
            "QfO mean is a project-defined unweighted six-metric secondary summary, not an official score.",
            "QfO metric_x means precision/recall or assessed relation count depending on the recorded axis, not total prediction coverage.",
            "Three Kingdoms scores only pairs among BUSCO reference genes, ignoring false positives involving other genes.",
            "OrthoBench non-primary-tool rows are from the prior scoring audit; their raw-output provenance is not yet consolidated here.",
            "Runtime budgets and memory accounting differ. Listed Three Kingdoms timings are historical evidence, not a matched efficiency experiment.",
            "FastOMA used supplied trees; tree-construction cost and reference-resource overlap require explicit accounting.",
            "Independent validation, controlled ablations, error strata, robustness, biological application, and release work remain open.",
        ],
    }


def markdown(report):
    lines = ["# Retained Benchmark Comparison", "", "Work in progress; not a frozen publication baseline.", "",
             "| Method | OrthoBench F1 (%) | QfO six-metric mean | Three Kingdoms F1 |",
             "| --- | ---: | ---: | ---: |"]
    for record in report["methods"]:
        qfo = record["qfo"]
        value = f"{qfo['mean']:.6f}" if qfo["status"] == "metrics_available" else "pending"
        lines.append(f"| {record['method']} | {record['orthobench']['f_score_percent']:.6f} | {value} | "
                     f"{record['three_kingdoms']['score']['f_score']:.6f} |")
    lines += ["", "## QfO Components", "", "| Method | " + " | ".join(METRICS) + " | Output |",
              "| --- | " + " | ".join(["---:"] * len(METRICS)) + " | --- |"]
    for record in report["methods"]:
        qfo = record["qfo"]
        scores = [f"{qfo['scores'][m]:.6f}" if qfo["status"] == "metrics_available" else "pending" for m in METRICS]
        lines.append("| " + record["method"] + " | " + " | ".join(scores) + " | " + qfo["output_semantics"] + " |")
    diagnostic = report["diagnostics"]["orthomcl_preclustering"]
    if "mean" in diagnostic:
        lines += ["", f"OrthoMCL pre-clustering diagnostic mean: {diagnostic['mean']:.6f}. "
                  "This is not the final-group assessment and is excluded from the main comparison."]
    lines += ["", "## Limitations", "", *["- " + note for note in report["limitations"]]]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = build_report(args.root)
    for path in (args.json, args.markdown):
        path.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    args.markdown.write_text(markdown(result))
    print(f"Wrote {len(result['methods'])} retained method rows; publication_ready=False")


if __name__ == "__main__":
    main()
