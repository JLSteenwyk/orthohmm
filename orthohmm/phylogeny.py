"""Configuration and species-tree validation for phylogeny-aware inference."""

from __future__ import annotations

from dataclasses import dataclass
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
    root_children = tree.seed_node.num_child_nodes()
    if explicitly_unrooted or root_children != 2:
        raise PhylogenyConfigurationError(
            "Species tree must be rooted and bifurcate at the root. For an "
            "unmarked rooted Newick tree, place the root on an edge so it has "
            "exactly two children; '[&U]' trees are not accepted."
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
