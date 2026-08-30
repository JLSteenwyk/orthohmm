"""Configuration and species-tree validation for phylogeny-aware inference."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import shutil
from typing import Sequence


FASTA_EXTENSIONS = (".fasta", ".prot", ".faa", ".fas", ".pep", ".fa")


class PhylogenyConfigurationError(ValueError):
    """Raised when a phylogeny option, dependency, or tree is invalid."""


@dataclass(frozen=True)
class PhylogenyConfig:
    mode: str = "off"
    species_tree_mode: str = "supplied"
    species_tree: str | None = None
    aligner: str = "mafft"
    tree_builder: str = "FastTree"


@dataclass(frozen=True)
class ValidatedSpeciesTree:
    """A parsed tree plus the canonical taxa used by OrthoHMM."""

    tree: object
    path: Path
    taxa: tuple[str, ...]


@dataclass(frozen=True)
class ReconciledNode:
    node_id: str
    parent_node_id: str | None
    event: str
    species_tree_node: str
    genes: tuple[str, ...]
    species: tuple[str, ...]


@dataclass(frozen=True)
class ReconciliationResult:
    family_id: str
    nodes: tuple[ReconciledNode, ...]
    root_groups: tuple[tuple[str, ...], ...]
    ortholog_pairs: tuple[tuple[str, str], ...]
    paralog_pairs: tuple[tuple[str, str], ...]
    duplication_count: int
    speciation_count: int
    annotated_newick: str


def canonical_species_name(filename: str) -> str:
    """Return the taxon label represented by an input proteome filename."""

    name = Path(filename).name
    lower_name = name.lower()
    for extension in FASTA_EXTENSIONS:
        if lower_name.endswith(extension):
            return name[: -len(extension)]
    return name


def canonical_species_names(proteome_files: Sequence[str]) -> tuple[str, ...]:
    taxa = tuple(canonical_species_name(filename) for filename in proteome_files)
    duplicates = sorted({taxon for taxon in taxa if taxa.count(taxon) > 1})
    if duplicates:
        raise PhylogenyConfigurationError(
            "Input proteome filenames do not map to unique taxon names: "
            + ", ".join(duplicates)
        )
    return taxa


def _load_dendropy():
    try:
        import dendropy
    except ImportError as exc:
        raise PhylogenyConfigurationError(
            "Phylogeny-aware inference requires the optional 'dendropy' "
            "dependency. Install OrthoHMM with 'orthohmm[phylogeny]'."
        ) from exc
    return dendropy


def parse_gene_tree(newick: str):
    """Parse one rooted gene tree while preserving exact gene identifiers."""

    dendropy = _load_dendropy()
    try:
        trees = dendropy.TreeList.get(
            data=newick,
            schema="newick",
            preserve_underscores=True,
            rooting="default-rooted",
        )
    except Exception as exc:
        raise PhylogenyConfigurationError(f"Invalid Newick gene tree: {exc}") from exc
    if len(trees) != 1:
        raise PhylogenyConfigurationError(
            f"Gene-tree input must contain exactly one tree; found {len(trees)}"
        )
    tree = trees[0]
    if tree.seed_node.num_child_nodes() < 2:
        raise PhylogenyConfigurationError(
            "Gene tree must contain at least two lineages at its root."
        )
    tree.is_rooted = True
    return tree


def _species_node_labels(species_tree) -> dict[object, str]:
    labels = {}
    internal_index = 0
    for node in species_tree.preorder_node_iter():
        if node.is_leaf():
            labels[node] = str(node.taxon.label)
        else:
            labels[node] = f"S{internal_index:04d}"
            internal_index += 1
    return labels


def reconcile_gene_tree(
    gene_tree,
    species_tree,
    gene_to_species: dict[str, str],
    family_id: str = "family",
) -> ReconciliationResult:
    """Reconcile a rooted gene tree by LCA mapping onto a rooted species tree."""

    species_leaves = {
        str(node.taxon.label): node for node in species_tree.leaf_node_iter()
    }
    species_labels = _species_node_labels(species_tree)
    species_depth = {
        node: len(list(node.ancestor_iter()))
        for node in species_tree.preorder_node_iter()
    }
    ancestor_paths = {
        label: (node,) + tuple(node.ancestor_iter())
        for label, node in species_leaves.items()
    }

    def species_lca(species: set[str]):
        unknown = sorted(species - species_leaves.keys())
        if unknown:
            raise PhylogenyConfigurationError(
                "Gene tree references species absent from the species tree: "
                + ", ".join(unknown)
            )
        common = set(ancestor_paths[next(iter(species))])
        for label in species:
            common.intersection_update(ancestor_paths[label])
        return max(common, key=species_depth.__getitem__)

    leaf_gene = {}
    observed_genes = []
    for leaf in gene_tree.leaf_node_iter():
        gene = str(leaf.taxon.label) if leaf.taxon is not None else ""
        if not gene:
            raise PhylogenyConfigurationError(
                "Every gene-tree leaf must have a non-empty gene identifier."
            )
        if gene not in gene_to_species:
            raise PhylogenyConfigurationError(
                f"Gene-tree leaf {gene!r} is absent from the candidate family."
            )
        leaf_gene[leaf] = gene
        observed_genes.append(gene)
    duplicate_genes = sorted(
        {gene for gene in observed_genes if observed_genes.count(gene) > 1}
    )
    if duplicate_genes:
        raise PhylogenyConfigurationError(
            "Gene-tree leaf identifiers are not unique: "
            + ", ".join(duplicate_genes)
        )

    descendant_genes = {}
    descendant_species = {}
    mapped_species_node = {}
    event_by_node = {}
    for node in gene_tree.postorder_node_iter():
        if node.is_leaf():
            gene = leaf_gene[node]
            genes = (gene,)
            species = {canonical_species_name(gene_to_species[gene])}
            event = "leaf"
        else:
            children = node.child_nodes()
            genes = tuple(
                gene
                for child in children
                for gene in descendant_genes[child]
            )
            child_species = [descendant_species[child] for child in children]
            species = set().union(*child_species)
            overlap = any(
                left & right
                for left, right in combinations(child_species, 2)
            )
            mapped = species_lca(species)
            child_maps_to_parent = any(
                mapped_species_node[child] is mapped for child in children
            )
            event = (
                "duplication" if overlap or child_maps_to_parent else "speciation"
            )
        descendant_genes[node] = tuple(sorted(genes))
        descendant_species[node] = species
        mapped_species_node[node] = species_lca(species)
        event_by_node[node] = event

    ortholog_pairs = set()
    paralog_pairs = set()
    for node in gene_tree.postorder_node_iter():
        if node.is_leaf():
            continue
        destination = (
            ortholog_pairs
            if event_by_node[node] == "speciation"
            else paralog_pairs
        )
        for left, right in combinations(node.child_nodes(), 2):
            for gene_a in descendant_genes[left]:
                for gene_b in descendant_genes[right]:
                    if gene_to_species[gene_a] == gene_to_species[gene_b]:
                        continue
                    destination.add(tuple(sorted((gene_a, gene_b))))

    species_root = species_tree.seed_node

    def root_lineages(node) -> list[tuple[str, ...]]:
        if (
            not node.is_leaf()
            and event_by_node[node] == "duplication"
            and mapped_species_node[node] is species_root
        ):
            groups = []
            for child in node.child_node_iter():
                groups.extend(root_lineages(child))
            return groups
        return [descendant_genes[node]]

    root_groups = tuple(
        sorted(root_lineages(gene_tree.seed_node), key=lambda group: group)
    )
    node_ids = {}
    internal_index = 0
    for node in gene_tree.postorder_node_iter():
        if node.is_leaf():
            node_ids[node] = leaf_gene[node]
        else:
            node_ids[node] = f"G{internal_index:05d}"
            internal_index += 1
    reconciled_nodes = []
    for node in gene_tree.postorder_node_iter():
        node_id = node_ids[node]
        if not node.is_leaf():
            event_code = "D" if event_by_node[node] == "duplication" else "S"
            node.label = (
                f"{node_id}|{event_code}@{species_labels[mapped_species_node[node]]}"
            )
        reconciled_nodes.append(
            ReconciledNode(
                node_id=node_id,
                parent_node_id=(
                    node_ids[node.parent_node]
                    if node.parent_node is not None
                    else None
                ),
                event=event_by_node[node],
                species_tree_node=species_labels[mapped_species_node[node]],
                genes=descendant_genes[node],
                species=tuple(sorted(descendant_species[node])),
            )
        )

    annotated_newick = gene_tree.as_string(
        schema="newick",
        suppress_rooting=False,
        suppress_internal_node_labels=False,
        preserve_spaces=True,
    ).strip()
    return ReconciliationResult(
        family_id=family_id,
        nodes=tuple(reconciled_nodes),
        root_groups=root_groups,
        ortholog_pairs=tuple(sorted(ortholog_pairs)),
        paralog_pairs=tuple(sorted(paralog_pairs)),
        duplication_count=sum(
            event == "duplication" for event in event_by_node.values()
        ),
        speciation_count=sum(
            event == "speciation" for event in event_by_node.values()
        ),
        annotated_newick=annotated_newick,
    )


def parse_species_tree(
    path: str | Path,
    proteome_files: Sequence[str],
) -> ValidatedSpeciesTree:
    """Parse and validate a rooted species tree against input proteomes."""

    tree_path = Path(path)
    if not tree_path.is_file():
        raise PhylogenyConfigurationError(
            f"Species tree does not exist or is not a file: {tree_path}"
        )

    try:
        newick = tree_path.read_text(encoding="utf-8").strip()
    except OSError as exc:
        raise PhylogenyConfigurationError(
            f"Could not read species tree {tree_path}: {exc}"
        ) from exc
    if not newick:
        raise PhylogenyConfigurationError(f"Species tree is empty: {tree_path}")

    dendropy = _load_dendropy()
    try:
        trees = dendropy.TreeList.get(
            data=newick,
            schema="newick",
            preserve_underscores=True,
            rooting="default-rooted",
        )
    except (
        dendropy.dataio.newickreader.NewickReader.NewickReaderDuplicateTaxonError
    ) as exc:
        raise PhylogenyConfigurationError(
            f"Species-tree taxon labels are not unique: {exc}"
        ) from exc
    except Exception as exc:
        raise PhylogenyConfigurationError(
            f"Invalid Newick species tree {tree_path}: {exc}"
        ) from exc
    if len(trees) != 1:
        raise PhylogenyConfigurationError(
            f"Species tree file must contain exactly one tree; found {len(trees)}"
        )

    tree = trees[0]
    explicitly_unrooted = newick.lstrip().upper().startswith("[&U]")
    explicitly_rooted = newick.lstrip().upper().startswith("[&R]")
    root_children = tree.seed_node.num_child_nodes()
    if explicitly_unrooted or root_children < 2 or (
        not explicitly_rooted and root_children != 2
    ):
        raise PhylogenyConfigurationError(
            "Species tree must be rooted. Unmarked Newick is accepted only "
            "when the root has exactly two children; explicitly rooted '[&R]' "
            "polytomies are accepted, and '[&U]' trees are not."
        )
    tree.is_rooted = True

    tree_labels = []
    for leaf in tree.leaf_node_iter():
        label = leaf.taxon.label if leaf.taxon is not None else None
        if not label:
            raise PhylogenyConfigurationError(
                "Every species-tree leaf must have a non-empty taxon label."
            )
        canonical_label = canonical_species_name(label)
        leaf.taxon.label = canonical_label
        tree_labels.append(canonical_label)

    duplicates = sorted(
        {label for label in tree_labels if tree_labels.count(label) > 1}
    )
    if duplicates:
        raise PhylogenyConfigurationError(
            "Species-tree taxon labels are not unique: " + ", ".join(duplicates)
        )

    expected = set(canonical_species_names(proteome_files))
    observed = set(tree_labels)
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    if missing or extra:
        details = []
        if missing:
            details.append("missing from tree: " + ", ".join(missing))
        if extra:
            details.append("not present in input proteomes: " + ", ".join(extra))
        raise PhylogenyConfigurationError(
            "Species-tree taxa do not match input proteomes ("
            + "; ".join(details)
            + ")."
        )

    return ValidatedSpeciesTree(
        tree=tree,
        path=tree_path.resolve(),
        taxa=tuple(sorted(observed)),
    )


def resolve_executable(command: str, role: str) -> str:
    resolved = shutil.which(command)
    if resolved is None:
        raise PhylogenyConfigurationError(
            f"Phylogeny-aware inference could not find the {role} executable "
            f"{command!r}. Provide its command or path explicitly."
        )
    return str(Path(resolved).resolve())


def validate_phylogeny_config(
    config: PhylogenyConfig,
    proteome_files: Sequence[str],
) -> PhylogenyConfig:
    """Validate enabled phylogeny settings without affecting off-mode runs."""

    if config.mode == "off":
        return config
    if config.mode != "reconcile":
        raise PhylogenyConfigurationError(
            f"Unsupported phylogeny mode: {config.mode!r}"
        )
    if config.species_tree_mode not in {"supplied", "infer"}:
        raise PhylogenyConfigurationError(
            f"Unsupported species-tree mode: {config.species_tree_mode!r}"
        )

    canonical_species_names(proteome_files)
    species_tree = config.species_tree
    if config.species_tree_mode == "supplied":
        if not species_tree:
            raise PhylogenyConfigurationError(
                "--species_tree is required when --species_tree_mode supplied."
            )
        parse_species_tree(species_tree, proteome_files)
    elif species_tree:
        raise PhylogenyConfigurationError(
            "--species_tree cannot be used with --species_tree_mode infer."
        )

    return PhylogenyConfig(
        mode=config.mode,
        species_tree_mode=config.species_tree_mode,
        species_tree=species_tree,
        aligner=resolve_executable(config.aligner, "aligner"),
        tree_builder=resolve_executable(config.tree_builder, "tree builder"),
    )
