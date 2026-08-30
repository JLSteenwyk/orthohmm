"""Execution layer for optional phylogeny-aware orthogroup inference."""

from __future__ import annotations

from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time
from typing import Sequence
import statistics

from .helpers import get_all_fasta_entries
from .phylogeny import (
    PhylogenyConfig,
    PhylogenyConfigurationError,
    ReconciliationResult,
    canonical_species_name,
    parse_gene_tree,
    parse_species_tree,
    reconcile_gene_tree,
    validate_phylogeny_config,
    _load_dendropy,
)


@dataclass(frozen=True)
class FamilyOutcome:
    family_id: str
    genes: tuple[str, ...]
    root_groups: tuple[tuple[str, ...], ...]
    ortholog_pairs: tuple[tuple[str, str], ...]
    reconciliation: ReconciliationResult | None
    checkpoint_hit: bool


@dataclass(frozen=True)
class PhylogenyStageResult:
    candidate_families: int
    reconciled_families: int
    bypassed_families: int
    checkpoint_hits: int
    root_hogs: int
    ortholog_pairs: int
    duplications: int
    speciations: int
    species_tree_families: int
    species_tree_checkpoint_hit: bool
    output_directory: str


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_clusters(path: Path) -> list[tuple[str, ...]]:
    clusters = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            genes = tuple(sorted(line.split()))
            if genes:
                clusters.append(genes)
    return clusters


def _load_sequence_data(
    fasta_directory: str,
    files: Sequence[str],
) -> tuple[dict[str, str], dict[str, str]]:
    nested = get_all_fasta_entries(fasta_directory, list(files))
    sequences = {}
    gene_species = {}
    for filename in files:
        species = canonical_species_name(filename)
        for gene, sequence in nested[filename].items():
            if gene in sequences:
                raise PhylogenyConfigurationError(
                    f"Gene identifier {gene!r} occurs in more than one proteome. "
                    "Phylogeny-aware inference requires globally unique identifiers."
                )
            sequences[gene] = sequence
            gene_species[gene] = species
    return sequences, gene_species


def family_requires_reconciliation(
    genes: Sequence[str],
    gene_species: dict[str, str],
) -> bool:
    species_counts = Counter(gene_species[gene] for gene in genes)
    return (
        len(genes) >= 3
        and len(species_counts) >= 2
        and max(species_counts.values(), default=0) > 1
    )


def select_species_tree_families(
    family_records: Sequence[tuple[str, tuple[str, ...]]],
    gene_species: dict[str, str],
    sequences: dict[str, str],
    expected_species: Sequence[str],
    max_families: int = 200,
) -> list[tuple[str, tuple[str, ...]]]:
    """Select stable, high-occupancy single-copy families for tree inference."""

    minimum_taxa = max(2, min(len(expected_species), math.ceil(
        len(expected_species) * 0.5
    )))
    eligible = []
    for family_id, genes in family_records:
        species = [gene_species[gene] for gene in genes]
        if len(species) != len(set(species)) or len(species) < minimum_taxa:
            continue
        median_length = statistics.median(len(sequences[gene]) for gene in genes)
        eligible.append((family_id, genes, len(species), median_length))
    eligible.sort(key=lambda row: (-row[2], -row[3], row[0], row[1]))
    selected = eligible[:max_families]

    represented = {
        gene_species[gene] for _, genes, _, _ in selected for gene in genes
    }
    missing = set(expected_species) - represented
    if missing:
        selected_ids = {family_id for family_id, _, _, _ in selected}
        for record in eligible[max_families:]:
            family_id, genes, _, _ = record
            family_species = {gene_species[gene] for gene in genes}
            if family_id not in selected_ids and family_species & missing:
                selected.append(record)
                selected_ids.add(family_id)
                represented.update(family_species)
                missing = set(expected_species) - represented
                if not missing:
                    break
    if missing:
        raise PhylogenyConfigurationError(
            "Cannot infer a species tree because no selected single-copy "
            "families represent these taxa: " + ", ".join(sorted(missing))
        )
    if not selected:
        raise PhylogenyConfigurationError(
            "Cannot infer a species tree: no multi-species single-copy "
            "candidate families passed the occupancy filter."
        )
    return [(family_id, genes) for family_id, genes, _, _ in selected]


def _read_fasta(path: Path) -> dict[str, str]:
    records = {}
    header = None
    parts = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(parts)
                header = line[1:].split()[0]
                parts = []
            elif header is not None:
                parts.append(line)
    if header is not None:
        records[header] = "".join(parts)
    return records


def _infer_species_tree(
    family_records: Sequence[tuple[str, tuple[str, ...]]],
    sequences: dict[str, str],
    gene_species: dict[str, str],
    files: Sequence[str],
    config: PhylogenyConfig,
    phylogeny_directory: Path,
    cpu: int,
):
    expected_species = tuple(sorted(canonical_species_name(name) for name in files))
    selected = select_species_tree_families(
        family_records,
        gene_species,
        sequences,
        expected_species,
    )
    species_directory = phylogeny_directory / "species_tree_inference"
    rooted_path = phylogeny_directory / "species_tree.rooted.nwk"
    checkpoint_path = species_directory / "checkpoint.json"
    payload = {
        "schema_version": 1,
        "families": [
            [
                family_id,
                [[gene, gene_species[gene], sequences[gene]] for gene in genes],
            ]
            for family_id, genes in selected
        ],
        "taxa": list(expected_species),
        "aligner": config.aligner,
        "tree_builder": config.tree_builder,
        "alignment_mode": "auto_single_thread",
        "tree_model": "LG",
        "rooting": "midpoint",
    }
    input_hash = _sha256_text(
        json.dumps(payload, sort_keys=True, separators=(",", ":"))
    )
    checkpoint_hit = False
    if checkpoint_path.is_file() and rooted_path.is_file():
        try:
            checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
            checkpoint_hit = (
                checkpoint.get("status") == "complete"
                and checkpoint.get("input_sha256") == input_hash
                and checkpoint.get("species_tree_sha256") == _sha256_file(rooted_path)
            )
        except (OSError, ValueError):
            checkpoint_hit = False
    if checkpoint_hit:
        return parse_species_tree(rooted_path, files), {
            "families": len(selected),
            "checkpoint_hit": True,
            "selected_family_ids": [family_id for family_id, _ in selected],
        }

    def align_family(record):
        family_id, genes = record
        input_path = species_directory / "candidate_fastas" / f"{family_id}.faa"
        alignment_path = species_directory / "alignments" / f"{family_id}.faa"
        token_to_gene = _write_token_fasta(input_path, genes, sequences)
        _run_to_file(
            _aligner_command(config.aligner, input_path),
            alignment_path,
            "species-tree family alignment",
        )
        aligned_tokens = _read_fasta(alignment_path)
        if set(aligned_tokens) != set(token_to_gene):
            raise PhylogenyConfigurationError(
                f"Aligner changed or omitted identifiers for {family_id}."
            )
        widths = {len(sequence) for sequence in aligned_tokens.values()}
        if len(widths) != 1:
            raise PhylogenyConfigurationError(
                f"Aligner returned unequal sequence lengths for {family_id}."
            )
        by_species = {
            gene_species[token_to_gene[token]]: sequence
            for token, sequence in aligned_tokens.items()
        }
        return family_id, next(iter(widths)), by_species

    with ThreadPoolExecutor(max_workers=max(1, min(int(cpu), len(selected)))) as pool:
        alignments = list(pool.map(align_family, selected))

    concatenated = {species: [] for species in expected_species}
    for _, width, by_species in alignments:
        for species in expected_species:
            concatenated[species].append(by_species.get(species, "-" * width))
    species_token_map = {
        f"s{index:08d}": species
        for index, species in enumerate(expected_species)
    }
    supermatrix_lines = []
    for token, species in species_token_map.items():
        sequence = "".join(concatenated[species])
        if not sequence.strip("-"):
            raise PhylogenyConfigurationError(
                f"Species-tree supermatrix contains no residues for {species}."
            )
        supermatrix_lines.extend((f">{token}", sequence))
    supermatrix_path = species_directory / "species_tree_alignment.faa"
    _atomic_write(supermatrix_path, "\n".join(supermatrix_lines) + "\n")
    raw_tree_path = species_directory / "species_tree.raw.nwk"
    _run_to_file(
        _tree_command(config.tree_builder, supermatrix_path),
        raw_tree_path,
        "species-tree inference",
    )
    tree = _root_external_gene_tree(
        raw_tree_path.read_text(encoding="utf-8"),
        species_token_map,
    )
    rooted_newick = tree.as_string(
        schema="newick",
        suppress_rooting=False,
        suppress_internal_node_labels=True,
        preserve_spaces=True,
    ).strip()
    _atomic_write(rooted_path, rooted_newick + "\n")
    checkpoint = {
        "schema_version": 1,
        "status": "complete",
        "input_sha256": input_hash,
        "species_tree_sha256": _sha256_file(rooted_path),
        "selected_family_ids": [family_id for family_id, _ in selected],
        "supermatrix_sha256": _sha256_file(supermatrix_path),
    }
    _atomic_write(
        checkpoint_path,
        json.dumps(checkpoint, indent=2, sort_keys=True) + "\n",
    )
    return parse_species_tree(rooted_path, files), {
        "families": len(selected),
        "checkpoint_hit": False,
        "selected_family_ids": [family_id for family_id, _ in selected],
        "supermatrix_columns": sum(
            len(part) for part in next(iter(concatenated.values()))
        ),
    }


def _write_token_fasta(
    path: Path,
    genes: Sequence[str],
    sequences: dict[str, str],
) -> dict[str, str]:
    token_to_gene = {
        f"g{index:08d}": gene for index, gene in enumerate(sorted(genes))
    }
    lines = []
    for token, gene in token_to_gene.items():
        lines.append(f">{token}")
        sequence = sequences[gene]
        lines.extend(sequence[i:i + 80] for i in range(0, len(sequence), 80))
    _atomic_write(path, "\n".join(lines) + "\n")
    return token_to_gene


def _run_to_file(command: list[str], output_path: Path, role: str) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    environment = os.environ.copy()
    environment.update({
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    })
    with temporary.open("w", encoding="utf-8") as output:
        completed = subprocess.run(
            command,
            stdout=output,
            stderr=subprocess.PIPE,
            text=True,
            env=environment,
            check=False,
        )
    if completed.returncode != 0:
        temporary.unlink(missing_ok=True)
        stderr = completed.stderr.strip()
        raise RuntimeError(
            f"{role} failed with exit code {completed.returncode}: "
            f"{' '.join(command)}\n{stderr}"
        )
    if not temporary.stat().st_size:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"{role} produced an empty output: {' '.join(command)}")
    os.replace(temporary, output_path)


def _aligner_command(executable: str, input_path: Path) -> list[str]:
    if "mafft" in Path(executable).name.lower():
        return [executable, "--thread", "1", "--auto", str(input_path)]
    return [executable, str(input_path)]


def _tree_command(executable: str, alignment_path: Path) -> list[str]:
    if "fasttree" in Path(executable).name.lower():
        return [executable, "-lg", str(alignment_path)]
    return [executable, str(alignment_path)]


def _root_external_gene_tree(newick: str, token_to_gene: dict[str, str]):
    explicitly_rooted = newick.lstrip().upper().startswith("[&R]")
    if explicitly_rooted:
        tree = parse_gene_tree(newick)
    else:
        dendropy = _load_dendropy()
        try:
            tree = dendropy.Tree.get(
                data=newick,
                schema="newick",
                preserve_underscores=True,
                rooting="default-unrooted",
            )
            tree.reroot_at_midpoint(update_bipartitions=False)
            tree.is_rooted = True
        except Exception as exc:
            raise PhylogenyConfigurationError(
                f"Could not parse and midpoint-root inferred gene tree: {exc}"
            ) from exc
    observed_tokens = set()
    for leaf in tree.leaf_node_iter():
        token = str(leaf.taxon.label) if leaf.taxon is not None else ""
        if token not in token_to_gene:
            raise PhylogenyConfigurationError(
                f"Tree builder returned unknown sequence identifier {token!r}."
            )
        observed_tokens.add(token)
        leaf.taxon.label = token_to_gene[token]
    missing = sorted(token_to_gene.keys() - observed_tokens)
    if missing:
        raise PhylogenyConfigurationError(
            "Tree builder omitted candidate sequences: " + ", ".join(missing)
        )
    return tree


def _tool_version(executable: str) -> str:
    basename = Path(executable).name.lower()
    command = [executable, "--version"] if "mafft" in basename else [executable]
    try:
        completed = subprocess.run(
            command,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return f"version unavailable: {exc}"
    text = "\n".join(part.strip() for part in (completed.stdout, completed.stderr))
    lines = [line for line in text.splitlines() if line.strip()]
    return lines[0][:500] if lines else "version unavailable"


def _family_input_hash(
    family_id: str,
    genes: Sequence[str],
    sequences: dict[str, str],
    gene_species: dict[str, str],
    species_tree_sha256: str,
    config: PhylogenyConfig,
) -> str:
    payload = {
        "schema_version": 1,
        "family_id": family_id,
        "sequences": [
            [gene, gene_species[gene], sequences[gene]] for gene in sorted(genes)
        ],
        "species_tree_sha256": species_tree_sha256,
        "aligner": config.aligner,
        "tree_builder": config.tree_builder,
        "alignment_mode": "auto_single_thread",
        "tree_model": "LG",
        "rooting": "explicit_or_midpoint",
    }
    return _sha256_text(json.dumps(payload, sort_keys=True, separators=(",", ":")))


def _reconcile_family(
    family_id: str,
    genes: tuple[str, ...],
    sequences: dict[str, str],
    gene_species: dict[str, str],
    species_tree,
    species_tree_sha256: str,
    config: PhylogenyConfig,
    phylogeny_directory: Path,
) -> FamilyOutcome:
    candidate_path = phylogeny_directory / "candidate_fastas" / f"{family_id}.faa"
    alignment_path = phylogeny_directory / "alignments" / f"{family_id}.faa"
    raw_tree_path = phylogeny_directory / "gene_trees" / f"{family_id}.raw.nwk"
    rooted_tree_path = phylogeny_directory / "gene_trees" / f"{family_id}.rooted.nwk"
    annotated_tree_path = (
        phylogeny_directory / "gene_trees" / f"{family_id}.reconciled.nwk"
    )
    checkpoint_path = phylogeny_directory / "checkpoints" / f"{family_id}.json"
    input_hash = _family_input_hash(
        family_id,
        genes,
        sequences,
        gene_species,
        species_tree_sha256,
        config,
    )

    checkpoint_hit = False
    if checkpoint_path.is_file() and rooted_tree_path.is_file():
        try:
            checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
            checkpoint_hit = (
                checkpoint.get("status") == "complete"
                and checkpoint.get("input_sha256") == input_hash
                and checkpoint.get("rooted_tree_sha256")
                == _sha256_file(rooted_tree_path)
            )
        except (OSError, ValueError):
            checkpoint_hit = False

    if checkpoint_hit:
        tree = parse_gene_tree(rooted_tree_path.read_text(encoding="utf-8"))
    else:
        token_to_gene = _write_token_fasta(candidate_path, genes, sequences)
        _run_to_file(
            _aligner_command(config.aligner, candidate_path),
            alignment_path,
            "alignment",
        )
        _run_to_file(
            _tree_command(config.tree_builder, alignment_path),
            raw_tree_path,
            "gene-tree inference",
        )
        tree = _root_external_gene_tree(
            raw_tree_path.read_text(encoding="utf-8"),
            token_to_gene,
        )
        rooted_newick = tree.as_string(
            schema="newick",
            suppress_rooting=False,
            suppress_internal_node_labels=True,
            preserve_spaces=True,
        ).strip()
        _atomic_write(rooted_tree_path, rooted_newick + "\n")

    reconciliation = reconcile_gene_tree(
        tree,
        species_tree,
        {gene: gene_species[gene] for gene in genes},
        family_id=family_id,
    )
    _atomic_write(annotated_tree_path, reconciliation.annotated_newick + "\n")
    checkpoint = {
        "schema_version": 1,
        "status": "complete",
        "family_id": family_id,
        "input_sha256": input_hash,
        "rooted_tree_sha256": _sha256_file(rooted_tree_path),
        "annotated_tree_sha256": _sha256_file(annotated_tree_path),
        "genes": list(genes),
    }
    _atomic_write(
        checkpoint_path,
        json.dumps(checkpoint, indent=2, sort_keys=True) + "\n",
    )
    return FamilyOutcome(
        family_id=family_id,
        genes=genes,
        root_groups=reconciliation.root_groups,
        ortholog_pairs=reconciliation.ortholog_pairs,
        reconciliation=reconciliation,
        checkpoint_hit=checkpoint_hit,
    )


def _bypass_family(
    family_id: str,
    genes: tuple[str, ...],
    gene_species: dict[str, str],
) -> FamilyOutcome:
    pairs = tuple(
        sorted(
            tuple(sorted((gene_a, gene_b)))
            for index, gene_a in enumerate(genes)
            for gene_b in genes[index + 1:]
            if gene_species[gene_a] != gene_species[gene_b]
        )
    )
    return FamilyOutcome(
        family_id=family_id,
        genes=genes,
        root_groups=(genes,),
        ortholog_pairs=pairs,
        reconciliation=None,
        checkpoint_hit=False,
    )


def _write_stage_outputs(
    outcomes: Sequence[FamilyOutcome],
    cluster_path: Path,
    phylogeny_directory: Path,
    gene_species: dict[str, str],
    manifest: dict,
) -> PhylogenyStageResult:
    root_rows = []
    replacement_clusters = []
    pairs = set()
    annotations = []
    hierarchical = []
    root_index = 0
    duplications = 0
    speciations = 0
    for outcome in outcomes:
        for group in outcome.root_groups:
            root_id = f"RootHOG{root_index:07d}"
            root_index += 1
            replacement_clusters.append(group)
            root_rows.append(
                "\t".join((root_id, outcome.family_id, ",".join(group)))
            )
        pairs.update(outcome.ortholog_pairs)
        result = outcome.reconciliation
        if result is None:
            hierarchical.append(
                "\t".join(
                    (f"{outcome.family_id}.root", "", "", outcome.family_id,
                     "unambiguous", ",".join(outcome.genes))
                )
            )
            continue
        duplications += result.duplication_count
        speciations += result.speciation_count
        for node in result.nodes:
            annotations.append(
                "\t".join(
                    (
                        outcome.family_id,
                        node.node_id,
                        node.parent_node_id or "",
                        node.event,
                        node.species_tree_node,
                        ",".join(node.species),
                        ",".join(node.genes),
                    )
                )
            )
            if node.event != "leaf":
                hierarchical.append(
                    "\t".join(
                        (
                            f"{outcome.family_id}.{node.node_id}",
                            (
                                f"{outcome.family_id}.{node.parent_node_id}"
                                if node.parent_node_id
                                else ""
                            ),
                            node.species_tree_node,
                            outcome.family_id,
                            node.event,
                            ",".join(node.genes),
                        )
                    )
                )

    cluster_text = "".join(
        " ".join(sorted(cluster)) + "\n" for cluster in replacement_clusters
    )
    _atomic_write(cluster_path, cluster_text)
    _atomic_write(
        phylogeny_directory / "orthohmm_root_hogs.tsv",
        "root_hog\tsource_family\tgenes\n" + "\n".join(root_rows) + "\n",
    )
    pair_rows = [
        "\t".join(
            (gene_a, gene_species[gene_a], gene_b, gene_species[gene_b])
        )
        for gene_a, gene_b in sorted(pairs)
    ]
    _atomic_write(
        phylogeny_directory / "orthohmm_pairwise_orthologs.tsv",
        "gene_a\tspecies_a\tgene_b\tspecies_b\n"
        + "\n".join(pair_rows)
        + "\n",
    )
    _atomic_write(
        phylogeny_directory / "orthohmm_reconciliation_nodes.tsv",
        "source_family\tnode_id\tparent_node_id\tevent\tspecies_tree_node\t"
        "species\tgenes\n"
        + "\n".join(annotations)
        + "\n",
    )
    _atomic_write(
        phylogeny_directory / "orthohmm_hierarchical_orthogroups.tsv",
        "hog_id\tparent_hog_id\tspecies_tree_node\tsource_family\tevent\tgenes\n"
        + "\n".join(hierarchical)
        + "\n",
    )

    reconciled = sum(outcome.reconciliation is not None for outcome in outcomes)
    checkpoint_hits = sum(outcome.checkpoint_hit for outcome in outcomes)
    result = PhylogenyStageResult(
        candidate_families=len(outcomes),
        reconciled_families=reconciled,
        bypassed_families=len(outcomes) - reconciled,
        checkpoint_hits=checkpoint_hits,
        root_hogs=len(replacement_clusters),
        ortholog_pairs=len(pairs),
        duplications=duplications,
        speciations=speciations,
        species_tree_families=int(
            manifest.get("species_tree_inference", {}).get("families", 0)
        ),
        species_tree_checkpoint_hit=bool(
            manifest.get("species_tree_inference", {}).get(
                "checkpoint_hit", False
            )
        ),
        output_directory=str(phylogeny_directory.resolve()),
    )
    summary = asdict(result)
    summary["schema_version"] = 1
    _atomic_write(
        phylogeny_directory / "reconciliation_summary.json",
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
    )
    manifest["results"] = summary
    _atomic_write(
        phylogeny_directory / "provenance_manifest.json",
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    return result


def run_phylogeny_stage(
    fasta_directory: str,
    output_directory: str,
    files: Sequence[str],
    config: PhylogenyConfig,
    cpu: int,
) -> PhylogenyStageResult:
    """Run supplied-tree reconciliation and replace clusters atomically."""

    config = validate_phylogeny_config(config, files)
    if config.mode != "reconcile":
        raise PhylogenyConfigurationError(
            "run_phylogeny_stage requires --phylogeny reconcile."
        )
    cluster_path = (
        Path(output_directory)
        / "orthohmm_working_res"
        / "orthohmm_edges_clustered.txt"
    )
    clusters = _read_clusters(cluster_path)
    sequences, gene_species = _load_sequence_data(fasta_directory, files)
    clustered_genes = {gene for cluster in clusters for gene in cluster}
    missing_sequences = sorted(clustered_genes - sequences.keys())
    if missing_sequences:
        raise PhylogenyConfigurationError(
            "Candidate clusters reference genes absent from the input proteomes: "
            + ", ".join(missing_sequences[:20])
        )

    phylogeny_directory = Path(output_directory) / "orthohmm_phylogeny"
    phylogeny_directory.mkdir(parents=True, exist_ok=True)
    family_records = [
        (f"Family{index:07d}", cluster) for index, cluster in enumerate(clusters)
    ]
    species_tree_inference = {"families": 0, "checkpoint_hit": False}
    if config.species_tree_mode == "supplied":
        validated_tree = parse_species_tree(str(config.species_tree), files)
        normalized_species_tree = validated_tree.tree.as_string(
            schema="newick",
            suppress_rooting=False,
            preserve_spaces=True,
        ).strip()
        species_tree_path = phylogeny_directory / "species_tree.rooted.nwk"
        _atomic_write(species_tree_path, normalized_species_tree + "\n")
        species_tree_source = str(validated_tree.path)
    else:
        validated_tree, species_tree_inference = _infer_species_tree(
            family_records,
            sequences,
            gene_species,
            files,
            config,
            phylogeny_directory,
            cpu,
        )
        species_tree_path = validated_tree.path
        species_tree_source = "internally_inferred_from_orthohmm_single_copy_families"
    species_tree_sha256 = _sha256_file(species_tree_path)
    ambiguous = [
        record
        for record in family_records
        if family_requires_reconciliation(record[1], gene_species)
    ]
    ambiguous_ids = {family_id for family_id, _ in ambiguous}
    outcomes_by_id = {
        family_id: _bypass_family(family_id, genes, gene_species)
        for family_id, genes in family_records
        if family_id not in ambiguous_ids
    }

    def worker(record):
        family_id, genes = record
        return _reconcile_family(
            family_id,
            genes,
            sequences,
            gene_species,
            validated_tree.tree,
            species_tree_sha256,
            config,
            phylogeny_directory,
        )

    with ThreadPoolExecutor(max_workers=max(1, min(int(cpu), len(ambiguous) or 1))) as pool:
        reconciled_outcomes = list(pool.map(worker, ambiguous))
    outcomes_by_id.update(
        {outcome.family_id: outcome for outcome in reconciled_outcomes}
    )
    outcomes = [outcomes_by_id[family_id] for family_id, _ in family_records]

    manifest = {
        "schema_version": 1,
        "created_at_epoch_s": time.time(),
        "mode": config.mode,
        "species_tree_mode": config.species_tree_mode,
        "species_tree_source": species_tree_source,
        "species_tree_sha256": species_tree_sha256,
        "species_tree_taxa": list(validated_tree.taxa),
        "species_tree_inference": species_tree_inference,
        "input_cluster_sha256": _sha256_file(cluster_path),
        "input_proteomes": [
            {
                "filename": filename,
                "taxon": canonical_species_name(filename),
                "sha256": _sha256_file(Path(fasta_directory) / filename),
            }
            for filename in files
        ],
        "cpu_budget": int(cpu),
        "parallelism": "deterministic_family_order_external_tools_one_thread_each",
        "gene_tree_rooting": "explicit_root_or_midpoint",
        "tools": {
            "aligner": {
                "path": config.aligner,
                "version": _tool_version(config.aligner),
            },
            "tree_builder": {
                "path": config.tree_builder,
                "version": _tool_version(config.tree_builder),
            },
        },
    }
    return _write_stage_outputs(
        outcomes,
        cluster_path,
        phylogeny_directory,
        gene_species,
        manifest,
    )
