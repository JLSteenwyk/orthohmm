"""Check candidate/alignment residues and reconstruct the species supermatrix."""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import Bio
from Bio import SeqIO

from benchmark_tools.audit_phylogeny_structure import rows, tree_leaves
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def sequences(path):
    result = {}
    for row in SeqIO.parse(path, "fasta"):
        sequence = str(row.seq)
        if not row.id or row.id in result or not sequence:
            raise ValueError(f"Empty or duplicate FASTA entry: {path}")
        result[row.id] = sequence
    if not result:
        raise ValueError(f"Empty FASTA: {path}")
    return result


def alignment_content(candidate, alignment, genes, originals):
    expected = {f"g{i:08d}": originals[gene] for i, gene in enumerate(sorted(genes))}
    if sequences(candidate) != expected:
        raise ValueError(f"Candidate tokens/residues differ from input: {candidate}")
    aligned = sequences(alignment)
    if set(aligned) != set(expected) or len({len(value) for value in aligned.values()}) != 1:
        raise ValueError(f"Alignment tokens or widths differ: {alignment}")
    for token, sequence in aligned.items():
        # MAFFT may change letter case; no residue substitution/deletion is allowed.
        if sequence.replace("-", "").upper() != expected[token].upper():
            raise ValueError(f"Alignment changed residues: {alignment}, {token}")
    return aligned


def audit(directory, structural_report):
    evidence = json.loads(structural_report.read_text())
    if evidence["status"] != "phylogeny_structure_verified":
        raise ValueError("Require successful structural readback")
    prior = {item["path"]: item for item in evidence["checked_records"]}
    watched = [record(structural_report)]

    def watch(path, prior_required=False):
        item = record(path)
        if prior_required and prior.get(item["path"]) != item:
            raise ValueError("Structural evidence no longer matches")
        watched.append(item)
        return item

    manifest_path = directory / "provenance_manifest.json"
    watch(manifest_path, True)
    manifest = json.loads(manifest_path.read_text())
    if manifest["species_tree_mode"] != "infer":
        raise ValueError("This audit requires inferred species-tree evidence")
    originals, gene_species = {}, {}
    for item in manifest["input_proteomes"]:
        matches = [Path(path) for path, rec in prior.items()
                   if Path(path).name == item["filename"] and rec["sha256"] == item["sha256"]]
        if len(matches) != 1:
            raise ValueError("Ambiguous input provenance")
        watch(matches[0], True)
        parsed = sequences(matches[0])
        if originals.keys() & parsed.keys():
            raise ValueError("Duplicate input genes")
        originals.update(parsed)
        gene_species.update(dict.fromkeys(parsed, item["taxon"]))
    roots = directory / "orthohmm_root_hogs.tsv"
    watch(roots, True)
    families = defaultdict(set)
    for row in rows(roots, ["root_hog", "source_family", "genes"]):
        families[row["source_family"]].update(row["genes"].split(","))
    reconciled = {Path(path).stem for path in prior
                  if Path(path).parent == directory / "checkpoints"}
    species_dir = directory / "species_tree_inference"
    cp_path = species_dir / "checkpoint.json"
    watch(cp_path)
    checkpoint = json.loads(cp_path.read_text())
    selected = checkpoint["selected_family_ids"]
    if (checkpoint["status"] != "complete" or not selected or len(selected) != len(set(selected))
            or selected != manifest["species_tree_inference"]["selected_family_ids"]
            or checkpoint["species_tree_sha256"] != manifest["species_tree_sha256"]):
        raise ValueError("Species-family checkpoint mismatch")
    concatenated = {taxon: [] for taxon in sorted(manifest["species_tree_taxa"])}
    counts = {}
    for stage, base, family_ids in (("gene", directory, sorted(reconciled)),
                                     ("species", species_dir, selected)):
        expected_names = {f"{family}.faa" for family in family_ids}
        for subdir in ("candidate_fastas", "alignments"):
            if {p.name for p in (base / subdir).iterdir()} != expected_names:
                raise ValueError(f"Unexpected {stage} {subdir} inventory")
        for family in family_ids:
            if family not in families:
                raise ValueError("Unknown alignment family")
            candidate = base / "candidate_fastas" / f"{family}.faa"
            alignment = base / "alignments" / f"{family}.faa"
            watch(candidate)
            watch(alignment)
            genes = sorted(families[family])
            aligned = alignment_content(candidate, alignment, genes, originals)
            if stage == "species":
                taxa = [gene_species[gene] for gene in genes]
                if len(taxa) != len(set(taxa)):
                    raise ValueError("Species marker is not single-copy")
                by_species = {taxon: aligned[f"g{i:08d}"] for i, taxon in enumerate(taxa)}
                width = len(next(iter(aligned.values())))
                for taxon in concatenated:
                    concatenated[taxon].append(by_species.get(taxon, "-" * width))
        counts[stage] = len(family_ids)
    supermatrix_path = species_dir / "species_tree_alignment.faa"
    if watch(supermatrix_path)["sha256"] != checkpoint["supermatrix_sha256"]:
        raise ValueError("Species supermatrix hash mismatch")
    expected = {f"s{i:08d}": "".join(parts) for i, parts in enumerate(concatenated.values())}
    if sequences(supermatrix_path) != expected:
        raise ValueError("Species supermatrix differs from ordered marker concatenation")
    raw_tree = species_dir / "species_tree.raw.nwk"
    watch(raw_tree)
    tree_leaves(raw_tree, expected)
    for item in watched:
        check(item)
    return dict(status="phylogeny_sequence_content_verified", gene_alignments=counts["gene"],
        species_alignments=counts["species"], supermatrix_columns=len(next(iter(expected.values()))),
        checked_records=watched, reader=record(__file__), biopython_version=Bio.__version__,
        scientific_scores_admitted=False, reconciliation_semantics_recomputed=False,
        limitations=["Exact candidate residues; aligned residues compared ignoring case and removing only '-' gaps.",
                     "Validates identity and concatenation, not biological alignment quality or tree optimality.",
                     "Node reconciliation events, pair completeness and native-run admission remain separate."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = audit(args.directory.resolve(), args.structure.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
