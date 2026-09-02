#!/usr/bin/env python3
"""Resume phylogeny after a verified fresh production run failed at that stage."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.native_pairs_to_qfo import convert
from benchmark_tools.orthobench_stage_diagnostics import (
    file_provenance,
    sha256_file,
)
from benchmark_tools.qfo_filter_pairs import filter_pairs, load_mapping
from orthohmm.files import fetch_fasta_files
from orthohmm.metrics import PipelineMetrics
from orthohmm.phylogeny import PhylogenyConfig
from orthohmm.phylogeny_pipeline import run_phylogeny_stage


CANDIDATE_RELATIVE_PATH = "orthohmm_working_res/orthohmm_edges_clustered.txt"


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def _git_commit(repo_root: Path) -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def verify_candidate_checkpoint(failed: dict, candidate: Path) -> dict:
    """Verify a candidate checkpoint against the failed harness manifest."""

    if failed.get("status") != "failed":
        raise ValueError("initial metrics must describe a failed production run")
    records = failed.get("harness", {}).get("output_manifest", [])
    matching = [
        record
        for record in records
        if record.get("path") == CANDIDATE_RELATIVE_PATH
    ]
    if len(matching) != 1:
        raise ValueError("failed metrics do not identify one candidate checkpoint")
    expected = matching[0]
    actual = file_provenance(candidate)
    if (
        actual["bytes"] != expected.get("bytes")
        or actual["sha256"] != expected.get("sha256")
    ):
        raise ValueError("candidate checkpoint does not match failed-run manifest")
    return actual


def combine_run_metrics(
    failed: dict,
    resumed: dict,
    stage_counts: dict,
) -> dict:
    """Combine the initial fresh run and its phylogeny restart conservatively."""

    stages = dict(failed.get("stages", {}))
    if "phylogeny" in stages:
        stages["phylogeny_initial_failure"] = stages.pop("phylogeny")
    stages["phylogeny_restart"] = resumed["stages"]["phylogeny"]
    counts = dict(failed.get("counts", {}))
    counts.update(stage_counts)
    initial_wall = float(failed.get("wall_s", 0.0))
    restart_wall = float(resumed.get("wall_s", 0.0))
    peak_bytes = max(
        int(failed.get("peak_process_tree_rss_bytes", 0)),
        int(resumed.get("peak_process_tree_rss_bytes", 0)),
    )
    return {
        "status": "complete",
        "wall_s": round(initial_wall + restart_wall, 6),
        "peak_process_tree_rss_bytes": peak_bytes,
        "peak_process_tree_rss_gib": round(peak_bytes / (1024 ** 3), 6),
        "stages": stages,
        "counts": counts,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta-directory", required=True, type=Path)
    parser.add_argument("--production-output", required=True, type=Path)
    parser.add_argument("--failed-metrics", required=True, type=Path)
    parser.add_argument("--resume-metrics", required=True, type=Path)
    parser.add_argument("--aggregate-json", required=True, type=Path)
    parser.add_argument("--qfo-mapping", required=True, type=Path)
    parser.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--aligner", default="mafft")
    parser.add_argument("--tree-builder", default="FastTree")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    repo_root = Path(__file__).resolve().parent.parent
    fasta_directory = args.fasta_directory.resolve()
    production_output = args.production_output.resolve()
    failed_metrics = args.failed_metrics.resolve()
    resume_metrics = args.resume_metrics.resolve()
    aggregate_json = args.aggregate_json.resolve()
    qfo_mapping = args.qfo_mapping.resolve()
    artifact_directory = aggregate_json.parent
    candidate = production_output / CANDIDATE_RELATIVE_PATH

    if aggregate_json.exists():
        raise SystemExit(f"aggregate result already exists: {aggregate_json}")
    failed = json.loads(failed_metrics.read_text(encoding="utf-8"))
    preserved_candidate = artifact_directory / "candidate_clusters.fresh_checkpoint.txt"
    try:
        if preserved_candidate.exists():
            candidate_record = verify_candidate_checkpoint(
                failed, preserved_candidate
            )
            if (
                not candidate.exists()
                or sha256_file(candidate) != candidate_record["sha256"]
            ):
                candidate.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(preserved_candidate, candidate)
        else:
            candidate_record = verify_candidate_checkpoint(failed, candidate)
            artifact_directory.mkdir(parents=True, exist_ok=True)
            shutil.copy2(candidate, preserved_candidate)
    except ValueError as error:
        raise SystemExit(str(error)) from error

    files = fetch_fasta_files(str(fasta_directory))
    config = PhylogenyConfig(
        mode="reconcile",
        species_tree_mode="infer",
        aligner=args.aligner,
        tree_builder=args.tree_builder,
    )
    restart_commit = _git_commit(repo_root)
    with PipelineMetrics(str(resume_metrics)) as metrics:
        metrics.add_metadata(
            harness="benchmark_tools.resume_phylogeny_benchmark",
            restart_kind="failed_fresh_end_to_end_checkpoint",
            initial_metrics=file_provenance(failed_metrics),
            candidate_checkpoint=file_provenance(preserved_candidate),
            initial_commit=failed.get("harness", {}).get("git_commit"),
            restart_commit=restart_commit,
            cpu_budget=args.cpu,
        )
        with metrics.stage("phylogeny"):
            stage_result = run_phylogeny_stage(
                str(fasta_directory),
                str(production_output),
                files,
                config,
                args.cpu,
            )
        metrics.add_counts(**asdict(stage_result))

    resumed = json.loads(resume_metrics.read_text(encoding="utf-8"))
    phylogeny_directory = production_output / "orthohmm_phylogeny"
    root_hogs = phylogeny_directory / "orthohmm_root_hogs.tsv"
    native_pairs = phylogeny_directory / "orthohmm_pairwise_orthologs.tsv"
    pairs = artifact_directory / "pairs.tsv"
    qfo_pairs = artifact_directory / "pairs.qfo.tsv"
    emitted_pairs = convert(native_pairs, pairs)
    total_pairs, retained_pairs = filter_pairs(
        pairs,
        qfo_pairs,
        load_mapping(qfo_mapping),
    )
    if emitted_pairs != total_pairs:
        raise RuntimeError("pair conversion and filtering counts disagree")
    shutil.copy2(root_hogs, artifact_directory / "root_hogs.tsv")

    combined = combine_run_metrics(failed, resumed, asdict(stage_result))
    combined.update({
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "run_kind": "fresh_end_to_end_with_verified_phylogeny_restart",
        "restart_reason": failed.get("error"),
        "information_budget": "internally_inferred_species_tree",
        "provenance": {
            "initial_commit": failed.get("harness", {}).get("git_commit"),
            "restart_commit": restart_commit,
            "initial_failed_metrics": file_provenance(failed_metrics),
            "restart_metrics": file_provenance(resume_metrics),
            "candidate_checkpoint": file_provenance(preserved_candidate),
            "qfo_mapping": file_provenance(qfo_mapping),
            "source": file_provenance(Path(__file__)),
        },
        "pair_preparation": {
            "native_pairs_emitted": emitted_pairs,
            "qfo_pairs_retained": retained_pairs,
            "qfo_pairs_removed": total_pairs - retained_pairs,
        },
        "outputs": {
            "species_tree": file_provenance(
                phylogeny_directory / "species_tree.rooted.nwk"
            ),
            "root_hogs": file_provenance(root_hogs),
            "native_pairs": file_provenance(native_pairs),
            "pairs": file_provenance(pairs),
            "qfo_pairs": file_provenance(qfo_pairs),
            "summary": file_provenance(
                phylogeny_directory / "reconciliation_summary.json"
            ),
            "manifest": file_provenance(
                phylogeny_directory / "provenance_manifest.json"
            ),
        },
    })
    _atomic_json(aggregate_json, combined)
    print(aggregate_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
