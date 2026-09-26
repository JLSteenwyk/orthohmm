# Phylogeny Sequence Readback

`benchmark_tools.audit_phylogeny_sequences` independently reads FASTAs with
Biopython 1.86 and consumes the prior structural readback. Input proteomes,
manifest and final source-family membership must match that evidence by hash.
The validator reconstructs token-to-gene associations from sorted family
members, checks exact candidate sequences, rejects missing/duplicate tokens
and unequal alignment widths, and compares ungapped aligned residues to the
input sequences. Only case differences and `-` alignment gaps are normalized;
residue substitutions or deletions are not accepted.

Both gene-tree and species-tree candidate/alignment directories must have
exactly the expected family files. Species-tree markers must be single-copy.
The recorded ordered marker list is used to reconstruct every species row,
including all-gap padding for absent taxa. The result must match the saved
supermatrix exactly; its checkpoint hash and raw species-tree token leaves
are also checked. All 17,780 full-dataset evidence records are rehashed at the
end of the readback.

## Results

- Installed synthetic fixture: passed; [receipt](installed_phylogeny_sequences_fixture_20260926.json).
- Historical OrthoBench p1_c1_r1: all 8,681 gene-family alignments and 200
  species-tree marker alignments pass. The 181,089-column supermatrix is
  reconstructed exactly. [Summary](installed_ob_historical_sequences_20260926.json).
- All 66 related reader and runner tests pass. Adversarial cases include token
  changes, substitutions, swapped candidates, unequal widths, duplicate FASTA
  entries, extra alignment files, duplicate markers and a hash-consistent but
  incorrectly concatenated supermatrix.

The full local historical report contains 17,780 file records (5.5 MB), SHA-256
`ef257beffe65e0047cdfe1bfdfb624b4c9e3fa096a8f51103fd72c062dd1fe50`.
It is not committed as a large generated manifest. Regenerate after the
[structural readback](PHYLOGENY_STRUCTURE_READBACK_20260926.md):

```bash
python -m benchmark_tools.audit_phylogeny_sequences \
  --directory benchmarks/results/publication_ob_factorial_v1/cells/p1_c1_r1/orthohmm_phylogeny \
  --structure benchmarks/work/installed_ob_historical_structure_20260926.json \
  --output benchmarks/work/installed_ob_historical_sequences_new.json
```

## Limits

Preserving residues does not establish biological alignment quality, tree
optimality, reconciliation correctness or completeness of emitted pairs. This
does not admit the still-running installed job 22179 or replace historical
scores. The frozen scientific method and immutable executor are unchanged.
Next validate node events and pair completeness against the saved trees and
membership policy, then apply all readback stages to the fresh run after native
completion. Other publication requirements remain open.
