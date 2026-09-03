#!/usr/bin/env python3
"""Replay the phylogeny stage from a frozen candidate-cluster checkpoint.

This development harness does not read benchmark labels during inference. It
records the candidate input, source revision, resource use, reconciliation
counts, and an optional post-run OrthoBench score in one JSON artifact.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import shutil
import sys
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.orthobench_stage_diagnostics import (
    file_provenance,
    git_state,
    run_official_benchmark,
)
from orthohmm.files import fetch_fasta_files
from orthohmm.metrics import PipelineMetrics
from orthohmm.phylogeny import PhylogenyConfig
from orthohmm.phylogeny_pipeline import run_phylogeny_stage


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta-directory", required=True, type=Path)
    parser.add_argument("--candidate-clusters", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    parser.add_argument(
        "--species-tree-mode", choices=("supplied", "infer"), default="infer"
    )
    parser.add_argument("--species-tree", type=Path)
    parser.add_argument("--aligner", default="mafft")
    parser.add_argument("--tree-builder", default="FastTree")
    parser.add_argument(
        "--root-rule",
        choices=(
            "supported_children", "confidence", "species_overlap", "mapped_event"
        ),
        default="supported_children",
    )
    parser.add_argument(
        "--pair-rule",
        choices=(
            "lca", "positive_paralogy", "supported_paralogy",
            "sparse_overlap", "depth_two_overlap",
        ),
        default="lca",
    )
    parser.add_argument(
        "--membership-support",
        choices=("high", "medium"),
        default="high",
    )
    parser.add_argument(
        "--membership-min-profile-support",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "--species-tree-rooting",
        choices=("midpoint", "min_variance"),
        default="midpoint",
    )
    parser.add_argument("--official-benchmark", type=Path)
    parser.add_argument(
        "--checkpoint-source",
        type=Path,
        help=(
            "Previous replay output whose immutable alignment and raw-tree "
            "checkpoints should seed this run"
        ),
    )
    parser.add_argument(
        "--membership-constraints",
        type=Path,
        help=(
            "Candidate-merge trace whose source and target seed groups should "
            "be verified after reconciliation"
        ),
    )
    return parser


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def _link_or_copy(source: str, destination: str) -> str:
    try:
        os.link(source, destination)
    except OSError:
        shutil.copy2(source, destination)
    return destination


def seed_checkpoint_output(source: Path, output: Path) -> Path:
    """Seed a new replay without allowing it to mutate the source run."""

    source = source.resolve()
    source_phylogeny = (
        source / "orthohmm_phylogeny"
        if (source / "orthohmm_phylogeny").is_dir()
        else source
    )
    required = (
        source_phylogeny / "checkpoints",
        source_phylogeny / "gene_trees",
    )
    missing = [str(path) for path in required if not path.is_dir()]
    if missing:
        raise SystemExit(
            "checkpoint source is missing required directories: "
            + ", ".join(missing)
        )
    destination = output.resolve() / "orthohmm_phylogeny"
    shutil.copytree(
        source_phylogeny,
        destination,
        copy_function=_link_or_copy,
    )
    return source_phylogeny


def load_membership_constraints(path: Path) -> list[dict]:
    """Load and validate candidate membership constraints for replay."""

    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, list):
        raise SystemExit("membership constraints must be a JSON list")
    for index, constraint in enumerate(payload):
        if not isinstance(constraint, dict):
            raise SystemExit(
                f"membership constraint {index} must be a JSON object"
            )
        missing = {"source_genes", "target_genes"} - constraint.keys()
        if missing:
            raise SystemExit(
                f"membership constraint {index} is missing: "
                + ", ".join(sorted(missing))
            )
    return payload


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.species_tree_mode == "supplied" and args.species_tree is None:
        raise SystemExit("--species-tree is required in supplied mode")
    if args.species_tree_mode == "infer" and args.species_tree is not None:
        raise SystemExit("--species-tree cannot be used in infer mode")

    output_directory = args.output_directory.resolve()
    cluster_path = (
        output_directory
        / "orthohmm_working_res"
        / "orthohmm_edges_clustered.txt"
    )
    if output_directory.exists() and any(output_directory.iterdir()):
        raise SystemExit(f"output directory is not empty: {output_directory}")
    checkpoint_source = None
    if args.checkpoint_source is not None:
        if args.checkpoint_source.resolve() == output_directory:
            raise SystemExit("checkpoint source and output directory must differ")
        checkpoint_source = seed_checkpoint_output(
            args.checkpoint_source,
            output_directory,
        )
    membership_constraints = None
    if args.membership_constraints is not None:
        membership_constraints = load_membership_constraints(
            args.membership_constraints.resolve()
        )
    cluster_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.candidate_clusters.resolve(), cluster_path)

    fasta_directory = args.fasta_directory.resolve()
    files = fetch_fasta_files(str(fasta_directory))
    config = PhylogenyConfig(
        mode="reconcile",
        species_tree_mode=args.species_tree_mode,
        species_tree=(
            str(args.species_tree.resolve()) if args.species_tree is not None else None
        ),
        aligner=args.aligner,
        tree_builder=args.tree_builder,
        root_duplication_rule=args.root_rule,
        pair_orthology_rule=args.pair_rule,
        species_tree_rooting=args.species_tree_rooting,
        membership_support_confidence=args.membership_support,
        membership_min_profile_support=args.membership_min_profile_support,
    )
    result_json = args.json.resolve()
    with PipelineMetrics(str(result_json)) as metrics:
        metrics.add_metadata(
            harness="benchmark_tools.replay_phylogeny",
            species_tree_mode=args.species_tree_mode,
            cpu_budget=args.cpu,
            root_duplication_rule=args.root_rule,
            pair_orthology_rule=args.pair_rule,
            membership_support_confidence=args.membership_support,
            membership_min_profile_support=args.membership_min_profile_support,
            species_tree_rooting=args.species_tree_rooting,
            membership_constraints=(
                file_provenance(args.membership_constraints)
                if args.membership_constraints is not None
                else None
            ),
        )
        with metrics.stage("phylogeny"):
            stage_result = run_phylogeny_stage(
                str(fasta_directory),
                str(output_directory),
                files,
                config,
                args.cpu,
                membership_constraints=membership_constraints,
            )
        metrics.add_counts(**asdict(stage_result))

    root_groups = (
        output_directory / "orthohmm_phylogeny" / "orthohmm_root_hogs.tsv"
    )
    payload = json.loads(result_json.read_text(encoding="utf-8"))
    payload.update({
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "git": git_state(),
        "source": file_provenance(Path(__file__)),
        "input": {
            "candidate_clusters": file_provenance(args.candidate_clusters),
            "fasta_directory": str(fasta_directory),
            "files": files,
        },
        "parameters": {
            "species_tree_mode": args.species_tree_mode,
            "species_tree": (
                file_provenance(args.species_tree)
                if args.species_tree is not None
                else None
            ),
            "aligner": args.aligner,
            "tree_builder": args.tree_builder,
            "root_duplication_rule": args.root_rule,
            "pair_orthology_rule": args.pair_rule,
            "species_tree_rooting": args.species_tree_rooting,
            "cpu": args.cpu,
            "checkpoint_source": (
                str(checkpoint_source) if checkpoint_source is not None else None
            ),
            "membership_constraints": (
                file_provenance(args.membership_constraints)
                if args.membership_constraints is not None
                else None
            ),
        },
        "outputs": {
            "root_hogs": file_provenance(root_groups),
            "summary": file_provenance(
                output_directory
                / "orthohmm_phylogeny"
                / "reconciliation_summary.json"
            ),
            "manifest": file_provenance(
                output_directory
                / "orthohmm_phylogeny"
                / "provenance_manifest.json"
            ),
        },
    })
    if args.official_benchmark is not None:
        payload["official_orthobench"] = run_official_benchmark(
            args.official_benchmark.resolve(), root_groups
        )
    _atomic_json(result_json, payload)
    score = payload.get("official_orthobench", {})
    print(
        f"root HOGs: F={score.get('f_score', 'NA')} "
        f"P={score.get('precision', 'NA')} R={score.get('recall', 'NA')}"
    )
    print(result_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
