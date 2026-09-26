# Phylogeny Structural Readback

Added `benchmark_tools.audit_phylogeny_structure` as an independent output
reader, using DendroPy 5.0.8 in the audit environment. It does not call the
production reconciliation implementation or alter predictions.

## Checks

- Input proteome identities, content hashes and unique gene-to-species mapping.
- Species-tree hash and exact taxa, final group coverage and source-family union.
- Complete checkpoint families against final source-family membership.
- Raw, rooted and reconciled tree hashes, one parsed tree per file, exact token
  or gene leaves, complete gene-tree directory inventory and family counts.
- Native pair and confidence tables agree row for row, have strict schemas,
  unique sorted canonical pairs, valid cross-species gene identities, matching
  source families, recognized confidence labels and the declared pair count.
- All recorded files are hashed again after readback to detect changes.

## Evidence

| Readback | Genes | Species | Source families | Reconciled families | Parsed trees | Native pairs |
|---|---:|---:|---:|---:|---:|---:|
| Installed synthetic fixture | 16 | 4 | 3 | 1 | 4 | 36 |
| Historical full OrthoBench p1_c1_r1 | 251,378 | 12 | 54,445 | 8,681 | 26,044 | 966,439 |

Historical final root groups: 59,770. Pair confidence: 849,010 high and 117,429
medium. These are internal consistency counts, not calibrated confidence
probabilities or new benchmark scores.

[Fixture receipt](installed_phylogeny_structure_fixture_20260926.json) and
[full-dataset summary](installed_ob_historical_structure_20260926.json) retain
reader identity and limitations. The full local 10.7 MB historical readback
contains 34,742 file records and has SHA-256
`0300c9e1462943be3efa810d8f6110e75d336e54cdcb930a90a9f389d6679464`;
it is not committed as a large generated manifest. Regenerate with:

```bash
python -m benchmark_tools.audit_phylogeny_structure \
  --directory benchmarks/results/publication_ob_factorial_v1/cells/p1_c1_r1/orthohmm_phylogeny \
  --input benchmarks/work/publication_installed_orthobench_20260926/input \
  --output benchmarks/work/installed_ob_historical_structure_new.json
```

The first historical invocation supplied a nonexistent original input path and
failed before readback. The corrected invocation used the existing byte-verified
private FASTA copies and passed. No native inference was retried or modified.
All 52 related structural-reader, installed-reader and runner tests pass.

## Scope Still Open

This is not admission of running job 22179. Its terminal-success gate, frozen
provenance, exact input inventory and score readback remain separate. Native
pairs are not equated with all pairs within root groups. Alignment/candidate
sequence content, node event semantics and completeness of emitted ortholog
pairs are not established by these structural checks. Full biological
correctness and publication readiness are not claimed. Preserve the prior
baseline score readback and run this reader on the new installed output only
after native completion.
