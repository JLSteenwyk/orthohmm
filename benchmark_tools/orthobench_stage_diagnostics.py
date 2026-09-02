#!/usr/bin/env python3
"""Measure where OrthoBench reference relationships disappear.

This is an offline diagnostic. It reads inference checkpoints after a run and
never exposes reference labels to the OrthoHMM inference process.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import pickle
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Sequence


GENE_PATTERN = re.compile(
    r"WBGene00\d+\.1|ENSCAFP\d+|ENSCINP\d+|ENSDARP\d+|FBpp0\d+|"
    r"ENSGALP\d+|ENSP000\d+|ENSMODP\d+|ENSMUSP\d+|ENSPTRP\d+|"
    r"ENSRNOP\d+|ENSTNIP\d+"
)


@dataclass(frozen=True)
class StageSpec:
    label: str
    path: Path


class _DisjointSet:
    def __init__(self, size: int):
        self.parent = list(range(size))
        self.size = [1] * size

    def find(self, item: int) -> int:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != item:
            parent = self.parent[item]
            self.parent[item] = root
            item = parent
        return root

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.size[left_root] < self.size[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        self.size[left_root] += self.size[right_root]


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_refogs(
    directory: Path, names: Sequence[str] | None = None
) -> list[set[str]]:
    paths = (
        [directory / name for name in names]
        if names is not None
        else sorted(directory.glob("RefOG*.txt"))
    )
    if not paths:
        raise ValueError(f"No RefOG*.txt files found in {directory}")
    groups = []
    for path in paths:
        genes = {
            line.strip()
            for line in path.read_text().splitlines()
            if line.strip()
        }
        if genes:
            groups.append(genes)
    if not groups:
        raise ValueError(f"No reference genes found in {directory}")
    return groups


def read_pair_file(path: Path) -> Iterator[tuple[str, str]]:
    """Read the first two whitespace-separated columns of a pair/ABC file."""
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected at least two columns")
            yield fields[0], fields[1]


def read_cached_hits(path: Path) -> Iterator[tuple[str, str]]:
    """Read pair keys from a local historical normalized-hit pickle."""
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    if not isinstance(payload, dict) or not isinstance(payload.get("all_hits"), dict):
        raise ValueError(f"{path} is not an OrthoHMM normalized-hit cache")
    for pair in payload["all_hits"]:
        if isinstance(pair, tuple) and len(pair) == 2:
            yield str(pair[0]), str(pair[1])


def read_clusters(path: Path) -> list[set[str]]:
    clusters = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            genes = set(GENE_PATTERN.findall(line))
            if not genes:
                fields = line.strip().split()
                if fields and fields[0].endswith(":"):
                    fields = fields[1:]
                genes = set(fields)
            if genes:
                clusters.append(genes)
    return clusters


def cluster_reference_pairs(
    clusters: Sequence[set[str]],
    reference_genes: set[str],
) -> Iterator[tuple[str, str]]:
    for cluster in clusters:
        retained = sorted(cluster & reference_genes)
        yield from itertools.combinations(retained, 2)


def relation_metrics(
    refogs: Sequence[set[str]],
    pairs: Iterable[tuple[str, str]],
) -> dict:
    gene_to_groups: dict[str, set[int]] = {}
    for group_index, group in enumerate(refogs):
        for gene in group:
            gene_to_groups.setdefault(gene, set()).add(group_index)
    reference_genes = set(gene_to_groups)
    local_indices = [
        {
            gene: index
            for index, gene in enumerate(sorted(group))
        }
        for group_index, group in enumerate(refogs)
    ]
    group_dsu = [_DisjointSet(len(group)) for group in refogs]
    direct_pairs: set[tuple[int, str, str]] = set()
    cross_pairs: set[tuple[str, str]] = set()
    observed_genes: set[str] = set()
    relation_records = 0

    for left, right in pairs:
        relation_records += 1
        if left == right or left not in reference_genes or right not in reference_genes:
            continue
        observed_genes.add(left)
        observed_genes.add(right)
        pair = tuple(sorted((left, right)))
        common_groups = gene_to_groups[left] & gene_to_groups[right]
        if not common_groups:
            cross_pairs.add(pair)
            continue
        for group_index in common_groups:
            group_pair = (group_index, *pair)
            if group_pair in direct_pairs:
                continue
            direct_pairs.add(group_pair)
            group_dsu[group_index].union(
                local_indices[group_index][left],
                local_indices[group_index][right],
            )

    total_pairs = 0
    transitive_pairs = 0
    macro_direct = []
    macro_transitive = []
    largest_component_fractions = []
    fully_connected = 0
    for group_index, group in enumerate(refogs):
        possible = len(group) * (len(group) - 1) // 2
        total_pairs += possible
        within_direct = sum(
            1
            for left, right in itertools.combinations(sorted(group), 2)
            if (group_index, left, right) in direct_pairs
        )
        component_sizes: dict[int, int] = {}
        dsu = group_dsu[group_index]
        for index in range(len(group)):
            root = dsu.find(index)
            component_sizes[root] = component_sizes.get(root, 0) + 1
        within_transitive = sum(
            size * (size - 1) // 2
            for size in component_sizes.values()
        )
        transitive_pairs += within_transitive
        if possible:
            macro_direct.append(within_direct / possible)
            macro_transitive.append(within_transitive / possible)
        largest = max(component_sizes.values(), default=0)
        largest_component_fractions.append(largest / len(group))
        if len(component_sizes) == 1:
            fully_connected += 1

    def mean(values: Sequence[float]) -> float:
        return sum(values) / len(values) if values else 0.0

    return {
        "relation_records": relation_records,
        "reference_genes": len(reference_genes),
        "reference_genes_observed": len(observed_genes),
        "reference_pairs": total_pairs,
        "direct_reference_pairs": len(direct_pairs),
        "direct_reference_pair_recall_micro": (
            len(direct_pairs) / total_pairs if total_pairs else 0.0
        ),
        "direct_reference_pair_recall_macro": mean(macro_direct),
        "transitive_reference_pairs": transitive_pairs,
        "transitive_reference_pair_recall_micro": (
            transitive_pairs / total_pairs if total_pairs else 0.0
        ),
        "transitive_reference_pair_recall_macro": mean(macro_transitive),
        "mean_largest_refog_component_fraction": mean(largest_component_fractions),
        "fully_connected_refogs": fully_connected,
        "refogs": len(refogs),
        "cross_refog_reference_pairs": len(cross_pairs),
    }


def parse_labeled_path(value: str) -> StageSpec:
    if "=" not in value:
        raise argparse.ArgumentTypeError("expected LABEL=PATH")
    label, raw_path = value.split("=", 1)
    if not label or not raw_path:
        raise argparse.ArgumentTypeError("expected non-empty LABEL=PATH")
    path = Path(raw_path)
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"file does not exist: {path}")
    return StageSpec(label, path)


def run_official_benchmark(script: Path, clusters: Path) -> dict:
    completed = subprocess.run(
        [sys.executable, str(script), str(clusters.resolve())],
        cwd=script.parent,
        check=True,
        capture_output=True,
        text=True,
    )
    metrics = {}
    patterns = {
        "f_score": r"([0-9.]+)% F-score",
        "precision": r"([0-9.]+)% Precision",
        "recall": r"([0-9.]+)% Recall",
        "exact_refogs": r"(\d+) Orthogroups exactly correct",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, completed.stdout)
        if match:
            metrics[key] = (
                int(match.group(1))
                if key == "exact_refogs"
                else float(match.group(1))
            )
    return metrics


def git_state() -> dict:
    def command(*args: str) -> str:
        completed = subprocess.run(
            ["git", *args], check=False, capture_output=True, text=True
        )
        return completed.stdout.strip() if completed.returncode == 0 else ""

    return {
        "commit": command("rev-parse", "HEAD") or None,
        "dirty": bool(command("status", "--porcelain")),
    }


def file_provenance(path: Path) -> dict:
    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refogs", type=Path, required=True)
    parser.add_argument(
        "--pairs", action="append", default=[], type=parse_labeled_path,
        metavar="LABEL=PATH", help="Two-column or ABC pair checkpoint",
    )
    parser.add_argument(
        "--hits-pickle", action="append", default=[], type=parse_labeled_path,
        metavar="LABEL=PATH", help="Trusted local historical all_hits pickle",
    )
    parser.add_argument(
        "--clusters", action="append", default=[], type=parse_labeled_path,
        metavar="LABEL=PATH", help="Orthogroup/cluster checkpoint",
    )
    parser.add_argument("--official-benchmark", type=Path)
    parser.add_argument("--partition-json", type=Path)
    parser.add_argument(
        "--partition", choices=("development", "validation", "all")
    )
    parser.add_argument("--json", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    labels = [spec.label for spec in args.pairs + args.hits_pickle + args.clusters]
    if len(labels) != len(set(labels)):
        raise SystemExit("Stage labels must be unique")
    if not labels:
        raise SystemExit("At least one stage checkpoint is required")

    if (args.partition_json is None) != (args.partition is None):
        raise SystemExit(
            "--partition-json and --partition must be specified together"
        )
    selected_names = None
    partition_provenance = None
    if args.partition_json is not None:
        partition_manifest = json.loads(args.partition_json.read_text())
        selected_names = (
            sorted(
                partition_manifest["development"]
                + partition_manifest["validation"]
            )
            if args.partition == "all"
            else partition_manifest[args.partition]
        )
        partition_provenance = file_provenance(args.partition_json)

    refogs = load_refogs(args.refogs, selected_names)
    reference_genes = set().union(*refogs)
    stages = []
    for spec in args.pairs:
        stages.append({
            "label": spec.label,
            "kind": "pairs",
            "input": file_provenance(spec.path),
            "metrics": relation_metrics(refogs, read_pair_file(spec.path)),
        })
    for spec in args.hits_pickle:
        stages.append({
            "label": spec.label,
            "kind": "cached_hits",
            "input": file_provenance(spec.path),
            "metrics": relation_metrics(refogs, read_cached_hits(spec.path)),
        })
    for spec in args.clusters:
        clusters = read_clusters(spec.path)
        record = {
            "label": spec.label,
            "kind": "clusters",
            "input": file_provenance(spec.path),
            "cluster_count": len(clusters),
            "metrics": relation_metrics(
                refogs, cluster_reference_pairs(clusters, reference_genes)
            ),
        }
        if args.official_benchmark:
            record["official_orthobench"] = run_official_benchmark(
                args.official_benchmark, spec.path
            )
        stages.append(record)

    output = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "command": [sys.executable, *sys.argv],
        "cwd": os.getcwd(),
        "git": git_state(),
        "source": file_provenance(Path(__file__)),
        "refogs": {
            "directory": str(args.refogs.resolve()),
            "partition": args.partition,
            "partition_manifest": partition_provenance,
            "groups": len(refogs),
            "genes": len(reference_genes),
            "manifest": [
                file_provenance(args.refogs / name)
                for name in (
                    selected_names
                    if selected_names is not None
                    else sorted(
                        path.name for path in args.refogs.glob("RefOG*.txt")
                    )
                )
            ],
        },
        "stages": stages,
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.json.with_suffix(args.json.suffix + ".tmp")
    temporary.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.json)

    print("| stage | direct pair recall | transitive pair recall | full RefOGs | cross-RefOG pairs |")
    print("| --- | ---: | ---: | ---: | ---: |")
    for stage in stages:
        metrics = stage["metrics"]
        print(
            f"| {stage['label']} | "
            f"{metrics['direct_reference_pair_recall_macro']:.3f} | "
            f"{metrics['transitive_reference_pair_recall_macro']:.3f} | "
            f"{metrics['fully_connected_refogs']}/{metrics['refogs']} | "
            f"{metrics['cross_refog_reference_pairs']} |"
        )
    print(args.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
