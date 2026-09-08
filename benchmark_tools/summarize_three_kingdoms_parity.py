#!/usr/bin/env python3
"""Build the machine-readable and Markdown Three Kingdoms parity reports."""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCORE_FIELDS = {
    "reference OGs": ("reference_orthogroups", int),
    "reference genes": ("reference_genes", int),
    "ref-genes in prediction": ("reference_genes_in_prediction", int),
    "predicted OGs": ("predicted_orthogroups", int),
    "TP gene pairs": ("true_positive_gene_pairs", int),
    "FP gene pairs": ("false_positive_gene_pairs", int),
    "FN gene pairs": ("false_negative_gene_pairs", int),
    "precision": ("precision", float),
    "recall": ("recall", float),
    "F-score": ("f_score", float),
}


@dataclass(frozen=True)
class MethodSpec:
    key: str
    tool: str
    version: str
    variant: str
    relative_result_dir: str
    phylogeny: bool
    cpus: int
    reused_inference: bool = False
    runtime_kind: str = "measured"
    memory_kind: str | None = "maximum RSS reported by GNU time"


METHODS = (
    MethodSpec(
        "orthohmm_high_sensitivity",
        "OrthoHMM",
        "0.5.0",
        "high sensitivity",
        "three_kingdoms/results/orthohmm_high_sensitivity",
        False,
        32,
        True,
        memory_kind="sampled peak process-tree RSS",
    ),
    MethodSpec(
        "orthohmm_phylogeny_satellite_v2",
        "OrthoHMM",
        "0.5.0",
        "high sensitivity + inferred phylogeny (satellite_v2)",
        "three_kingdoms/results/parity_20260907/orthohmm_phylogeny_satellite_v2",
        True,
        32,
        memory_kind="sampled peak process-tree RSS",
    ),
    MethodSpec(
        "orthofinder_3_1_5_full",
        "OrthoFinder",
        "3.1.5",
        "full pipeline; root HOGs",
        "three_kingdoms/results/parity_20260907/orthofinder_3_1_5_full",
        True,
        32,
    ),
    MethodSpec(
        "orthofinder_3_1_5_sequence_only",
        "OrthoFinder",
        "3.1.5",
        "sequence-only MCL checkpoint",
        "three_kingdoms/results/parity_20260907/orthofinder_3_1_5_sequence_only",
        False,
        32,
        runtime_kind="derived from matching full run's MCL checkpoint",
        memory_kind=None,
    ),
    MethodSpec(
        "proteinortho_6_3_6",
        "ProteinOrtho",
        "6.3.6",
        "default ProteinOrtho mode; DIAMOND",
        "three_kingdoms/results/proteinortho",
        False,
        32,
        True,
    ),
    MethodSpec(
        "sonicparanoid_2_0_9",
        "SonicParanoid",
        "2.0.9",
        "default mode; DIAMOND very-sensitive",
        "three_kingdoms/results/sonicparanoid",
        False,
        32,
        True,
    ),
    MethodSpec(
        "fastoma_0_3_5",
        "FastOMA",
        "0.3.5",
        "root HOGs; supplied species tree",
        "three_kingdoms/results/parity_20260907/fastoma_0_3_5",
        True,
        128,
    ),
    MethodSpec(
        "orthomcl_1_4",
        "OrthoMCL",
        "1.4",
        "defaults; legacy BLASTP; MCL inflation 1.5",
        "three_kingdoms/results/parity_20260907/orthomcl_1_4",
        False,
        32,
    ),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_score(path: Path) -> dict[str, int | float]:
    """Parse the stable text output from score_against_busco.py."""
    result: dict[str, int | float] = {}
    for line in path.read_text().splitlines():
        stripped = line.strip()
        for label, (key, converter) in SCORE_FIELDS.items():
            if not stripped.startswith(f"{label}"):
                continue
            raw = stripped.split(":", 1)[1].strip().split()[0].replace(",", "")
            result[key] = converter(raw)
            break
    missing = {key for key, _ in SCORE_FIELDS.values()} - result.keys()
    if missing:
        raise ValueError(f"Missing fields in {path}: {sorted(missing)}")
    reference_genes = int(result["reference_genes"])
    represented = int(result["reference_genes_in_prediction"])
    result["reference_gene_coverage"] = represented / reference_genes
    return result


def read_number(path: Path) -> int | float | None:
    if not path.exists():
        return None
    value = path.read_text().strip()
    if not value:
        return None
    number = float(value)
    return int(number) if number.is_integer() else number


def parse_run_metadata(path: Path) -> dict[str, int | str]:
    if not path.exists():
        return {}
    result: dict[str, int | str] = {}
    for line in path.read_text().splitlines():
        key, value = line.split("\t", 1)
        result[key] = int(value) if value.isdigit() else value
    return result


def read_runtime(result_dir: Path) -> int | float | None:
    for name in ("inference_wall_time_s.txt", "wall_time_s.txt"):
        value = read_number(result_dir / name)
        if value is not None:
            return value
    metrics = result_dir / "metrics.json"
    if metrics.exists():
        return json.loads(metrics.read_text()).get("wall_s")
    return None


def read_peak_rss_kib(result_dir: Path) -> int | float | None:
    for name in ("peak_rss_kib.txt", "peak_rss.txt"):
        value = read_number(result_dir / name)
        if value is not None:
            return value
    metrics = result_dir / "metrics.json"
    if metrics.exists():
        value = json.loads(metrics.read_text()).get("peak_process_tree_rss_bytes")
        if value is not None:
            return value / 1024
    return None


def collect_method(root: Path, spec: MethodSpec) -> dict[str, Any]:
    result_dir = root / spec.relative_result_dir
    score_path = result_dir / "score.txt"
    groups_path = result_dir / "orthogroups.txt"
    if not score_path.exists() or not groups_path.exists():
        missing = [str(path) for path in (score_path, groups_path) if not path.exists()]
        raise FileNotFoundError(", ".join(missing))

    metadata = parse_run_metadata(result_dir / "run_metadata.tsv")
    peak_rss_kib = read_peak_rss_kib(result_dir)
    return {
        "key": spec.key,
        "tool": spec.tool,
        "version": spec.version,
        "variant": spec.variant,
        "uses_phylogeny": spec.phylogeny,
        "result_dir": spec.relative_result_dir,
        "reused_completed_inference": spec.reused_inference,
        "score": parse_score(score_path),
        "performance": {
            "wall_s": read_runtime(result_dir),
            "runtime_kind": spec.runtime_kind,
            "cpus_requested": spec.cpus,
            "peak_rss_kib": peak_rss_kib,
            "peak_rss_gib": peak_rss_kib / 1024 / 1024 if peak_rss_kib else None,
            "memory_measurement": spec.memory_kind,
        },
        "provenance": {
            "slurm_job_id": metadata.get("slurm_job_id"),
            "source_commit": (
                (result_dir / "source_commit.txt").read_text().strip()
                if (result_dir / "source_commit.txt").exists()
                else None
            ),
            "orthogroups_sha256": sha256(groups_path),
            "score_sha256": sha256(score_path),
        },
    }


def format_duration(seconds: int | float | None) -> str:
    if seconds is None:
        return "n/a"
    seconds = round(seconds)
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:d}:{minutes:02d}:{seconds:02d}"


def render_markdown(report: dict[str, Any]) -> str:
    rows = sorted(report["methods"], key=lambda method: method["score"]["f_score"], reverse=True)
    lines = [
        "# Three Kingdoms Method Parity Benchmark",
        "",
        (
            "All methods were evaluated against the same 255 BUSCO reference "
            "orthogroups (2,035 reference genes) from the same 12-proteome, "
            "443,217-protein dataset."
        ),
        "",
        "| Method | Variant | Precision | Recall | F-score | Ref. coverage | Predicted OGs | Wall time | Peak RSS |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for method in rows:
        score = method["score"]
        perf = method["performance"]
        rss = perf["peak_rss_gib"]
        rss_text = f"{rss:.2f} GiB" if rss is not None else "n/a"
        lines.append(
            "| "
            f"{method['tool']} {method['version']} | {method['variant']} | "
            f"{score['precision']:.4f} | {score['recall']:.4f} | "
            f"{score['f_score']:.4f} | {score['reference_gene_coverage']:.1%} | "
            f"{score['predicted_orthogroups']:,} | "
            f"{format_duration(perf['wall_s'])} | {rss_text} |"
        )

    by_key = {method["key"]: method for method in report["methods"]}
    sensitive = by_key["orthohmm_high_sensitivity"]["score"]["f_score"]
    phylogeny = by_key["orthohmm_phylogeny_satellite_v2"]["score"]["f_score"]
    of_full = by_key["orthofinder_3_1_5_full"]["score"]["f_score"]
    of_sequence = by_key["orthofinder_3_1_5_sequence_only"]["score"]["f_score"]
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            (
                f"OrthoHMM phylogeny improves F-score by {phylogeny - sensitive:+.4f} "
                f"over high sensitivity alone and is {phylogeny - of_full:+.4f} versus "
                f"OrthoFinder's full root-HOG output. OrthoFinder's sequence-only "
                f"checkpoint remains higher by {of_sequence - phylogeny:.4f}."
            ),
            "",
            (
                "The benchmark is deliberately narrow: it scores gene-pair recovery only "
                "within conserved BUSCO families. Precision is often 1.0 because splitting "
                "a BUSCO family reduces recall while creating no cross-family false-positive "
                "pairs. This result should be interpreted alongside QfO and OrthoBench."
            ),
            "",
            "## Provenance Notes",
            "",
            "- OrthoFinder sequence-only time is derived from the matching full 3.1.5 run at its MCL checkpoint; it is not a separate timed invocation.",
            "- ProteinOrtho and SonicParanoid inference outputs were retained because they already used this exact biological input; they were rescored with the common evaluator for this report.",
            "- Root-HOG outputs are reported for phylogenetic pipelines, while flat orthogroups are reported for sequence-only methods.",
            "- External-tool GNU time RSS values may omit memory held by container descendants and are not directly comparable to OrthoHMM's sampled process-tree RSS.",
            "- Historical OrthoFinder 2.5.5 runs are excluded from the parity table; OrthoFinder 3.1.5 is the retained comparator.",
            "",
            "Machine-readable details, checksums, job IDs, and source revisions are in `three_kingdoms_parity_20260907.json`.",
            "",
        ]
    )
    return "\n".join(lines)


def build_report(root: Path) -> dict[str, Any]:
    reference = root / "three_kingdoms/busco/reference_orthogroups.txt"
    methods = [collect_method(root, spec) for spec in METHODS]
    return {
        "schema_version": 1,
        "benchmark": "Three Kingdoms BUSCO",
        "priority": "secondary",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "dataset": {
            "proteomes": 12,
            "proteins": 443217,
            "reference_orthogroups": 255,
            "reference_genes": 2035,
            "reference_sha256": sha256(reference),
        },
        "scoring": {
            "metric": "gene-pair micro-averaged precision, recall, and F-score",
            "script": "three_kingdoms/score_against_busco.py",
        },
        "methods": methods,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--json-output",
        type=Path,
        default=Path("benchmark_tools/results/three_kingdoms_parity_20260907.json"),
    )
    parser.add_argument(
        "--markdown-output",
        type=Path,
        default=Path("benchmark_tools/results/THREE_KINGDOMS_PARITY_20260907.md"),
    )
    args = parser.parse_args()

    root = args.root.resolve()
    report = build_report(root)
    json_output = args.json_output if args.json_output.is_absolute() else root / args.json_output
    markdown_output = (
        args.markdown_output
        if args.markdown_output.is_absolute()
        else root / args.markdown_output
    )
    json_output.parent.mkdir(parents=True, exist_ok=True)
    markdown_output.parent.mkdir(parents=True, exist_ok=True)
    json_output.write_text(json.dumps(report, indent=2) + "\n")
    markdown_output.write_text(render_markdown(report))
    print(f"Wrote {json_output}")
    print(f"Wrote {markdown_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
