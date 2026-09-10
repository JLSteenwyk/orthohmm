#!/usr/bin/env python3
"""Convert orthogroups.txt to QfO pairwise orthologs TSV.

Each line of orthogroups.txt is one OG: space-separated gene IDs.
Output (stdout): one row per ortholog pair, two columns separated by TAB.
A pair is emitted only when the two genes belong to *different* species.

Species assignment: each input proteome FASTA is one species; gene IDs are
mapped to species by reading the headers of every *.fasta in INPUT_DIR. This
matches QfO submission format expectations.

Usage:
    og_to_pairwise.py <orthogroups.txt> <input_proteomes_dir>
"""
import sys
from pathlib import Path
from itertools import combinations


def index_genes(input_dir):
    """Return dict gene_id -> species_filename.

    Indexes the FULL first-whitespace-separated token from each FASTA
    header (e.g. "sp|H2VFI5|SIAA_NEIMB"), since that is what every
    orthology tool emits in its orthogroup output. Also indexes the bare
    UniProt accession ("H2VFI5") as a fallback for tools that strip to
    the accession.
    """
    gene_to_species = {}
    for fa in sorted(input_dir.glob("*.fasta")):
        sp = fa.stem
        with open(fa) as f:
            for line in f:
                if line.startswith(">"):
                    tok = line[1:].split()[0]
                    gene_to_species[tok] = sp
                    if "|" in tok:
                        parts = tok.split("|")
                        if len(parts) >= 2:
                            gene_to_species.setdefault(parts[1], sp)
    return gene_to_species


def _strip_to_uniprot(token: str) -> str:
    """Reduce a FASTA header token to its bare UniProt accession.

    QfO's benchmark service rejects predictions that use UniProt-style
    composite identifiers like "sp|P12345|GENE_HUMAN" — it requires the
    bare accession ("P12345"). Returns the input unchanged when the token
    isn't pipe-separated.
    """
    if "|" in token:
        parts = token.split("|")
        if len(parts) >= 2 and parts[1]:
            return parts[1]
    return token


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: og_to_pairwise.py <orthogroups.txt> <input_proteomes_dir>")
    og_path = Path(sys.argv[1])
    input_dir = Path(sys.argv[2])
    if not og_path.exists():
        sys.exit(f"not found: {og_path}")
    if not input_dir.is_dir():
        sys.exit(f"not a directory: {input_dir}")

    gene_to_species = index_genes(input_dir)
    if not gene_to_species:
        sys.exit(f"no gene IDs found in {input_dir}/*.fasta")

    out = sys.stdout
    n_emitted = 0
    observed_genes = {}
    with open(og_path) as f:
        for line_number, line in enumerate(f, start=1):
            genes = line.split()
            if not genes:
                continue
            if len(genes) != len(set(genes)):
                raise ValueError(f"duplicate gene within {og_path}:{line_number}")
            for gene in genes:
                previous_line = observed_genes.get(gene)
                if previous_line is not None:
                    raise ValueError(
                        f"gene {gene!r} occurs in multiple groups in {og_path}: "
                        f"lines {previous_line} and {line_number}"
                    )
                if gene not in gene_to_species:
                    raise ValueError(
                        f"gene {gene!r} at {og_path}:{line_number} is absent from "
                        f"{input_dir}"
                    )
                observed_genes[gene] = line_number
            if len(genes) < 2:
                continue
            for a, b in combinations(genes, 2):
                sa = gene_to_species.get(a)
                sb = gene_to_species.get(b)
                if sa == sb:
                    continue
                a_acc = _strip_to_uniprot(a)
                b_acc = _strip_to_uniprot(b)
                # Canonicalize order so each pair appears once
                if a_acc > b_acc:
                    a_acc, b_acc = b_acc, a_acc
                out.write(f"{a_acc}\t{b_acc}\n")
                n_emitted += 1

    print(f"Emitted {n_emitted:,} pairs", file=sys.stderr)


if __name__ == "__main__":
    main()
