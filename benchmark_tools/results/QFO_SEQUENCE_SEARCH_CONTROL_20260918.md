# Corrected QfO Sequence-Search Control

## Scope

Extend the exploratory OrthoBench initial-search replacement to the frozen
corrected QfO2020_04 panel:984,137 proteins in78 species. This remains
development-exposed evidence, not independent confirmation. No defaults
change, and no QfO endpoint is used to select parameters.

Use the existing DIAMOND2.1.11 binary and unchanged search_command from
prepare_sequence_search_control.py:32 threads, very-sensitive, E-value1e-4,
BLOSUM62, gap11/1, composition1, masking1, all targets and one HSP. Query
all proteins against each species-specific database. Retain direction and
self hits; apply no identity/coverage cut. Preserve raw scores, bit scores,
E-values and lengths for audited normalization once by sqrt(query*target
length). Original OrthoBench manifest supplies binary/settings provenance
only; no original hit table or checkpoint is reused.

Primary comparison is P0/C0/R0 downstream graph processing versus the
corrected initial-HMM-search control, with the same frozen graph settings.
Retain all-hit and post-search top100 diagnostic variants as specified in
SEQUENCE_SEARCH_CONTROL_PROTOCOL_20260916.md. The post-search cap is not
the HMM prefilter cap. Equal E-values and identical graph parameters do
not establish equal sensitivity, score calibration or computational cost.
Do not claim that this alone solves the matched-sensitivity requirement.

## Prepared Inputs

Preparation implementation3cbe3a3 verifies pinned corrected staging and
direct inventory, every FASTA hash, exact gene/species scope and the existing
DIAMOND binary hash/version. It writes no search result and authorizes no
execution.24 focused preparation/conversion tests passed in0.26s.

Root: benchmarks/work/qfo_sequence_search_control_v1.

| Artifact | Bytes | SHA-256 |
| --- | ---: | --- |
| manifest.json | 205271 | fcef19d9c745c5806197aa219dec1c89bd5d90082a35f24dadafa3422f461436 |
| queries.fasta | 584039577 | 00e106d9e6f0eaf44cbaca0fbde6c63e2f587b848e6ca9cbaeb5d7a2268b5519 |
| gene_metadata.json | 80596465 | 3ce892e98e5254270729d44c187ac57fb706154385d590e3048a98c28f49b511 |

Post-preparation checks rehashed source, queries, metadata and binary;
all78 target directories are empty. Large derived files remain untracked.

## Remaining Gates

- Freeze a corrected-only runner with pre/post checks and full-panel
  success accounting; retain failure logs and forbid partial inference.
- Measure database and search phases separately. Shared-workstation cost
  is descriptive, not dedicated timing evidence. Do not load the active DGX.
- Validate hit ID ownership, exact lengths, unique ordered pairs, finite
  scores and target scope before numeric checkpoints. SQLite/memmap keeps
  conversion chunked but does not prove the downstream graph fits memory.
- Review observed output size and peak memory before graph replay; do not
  truncate all-hit data silently to meet resource limits. The historical
  OrthoBench all-hit checkpoint has100,099,147 hits, not a prediction of
  corrected QfO output size.
- Compare hit overlap/coverage and score distributions against the admitted
  corrected HMM checkpoint before accuracy interpretation, without fitting
  scale factors to benchmark outcomes. Independently validate and score
  both complete diagnostic variants, retaining negative findings.

No search has been launched and no new accuracy or timing result exists.
