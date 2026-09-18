# FastOMA Native Pair Conversion Audit

## Verified Result

The original-release FastOMA native `orthologs.tsv.gz` contains 15,320,615
pair rows. Strict parsing and accession ownership checks against all 78
frozen input FASTAs (976,504 accessions) pass. Canonical orientation of every
native row matches the retained `pairs.tsv` exactly, in order and with
multiplicity preserved. No malformed, foreign-accession, self, or
within-species pair was accepted. Both streams end at the same row.

This addresses the old converter's potential to silently discard malformed
rows: no such discrepancy was found in this retained output. Predictions
and scores were not modified.

Machine-readable evidence: `fastoma_pair_conversion_audit_20260918.json`,
SHA-256 `2cc3b75b0611dc1bf3a832b880ea0ad9cbe5ecd414fd0d77e8600b46339c2765`.
Source and all recorded input hashes were independently rechecked after the
audit. The source gzip checksum is
`c4cd1edaed5979a5a256630b19accee86b842528f0804b92af2903b097eba90f`;
retained TSV checksum is
`ffa1a85f12c0201b9a463fedf8b02c7ff08e9162a5760d8368227ca8063270cf`.

## Scope and Tests

This is not a corrected-release result or a validation of biological
orthology. Global pair uniqueness and original workflow completeness are
not asserted. The streaming converter preserves duplicate rows explicitly;
any downstream uniqueness requirement needs its own check.

The converter writes a stream and may have emitted a partial prefix before
rejecting a late malformed row. Callers must check its exit status and use a
temporary output before publishing a completed conversion.

Validation: 19 tests passed in
`test_fastoma_to_pairwise.py` and `test_audit_fastoma_pair_conversion.py`,
covering gzip/plain text, orientation, order, repeated rows, malformed
records, accession ownership, empty inputs/outputs, truncated gzip, and
missing/extra/changed retained rows.

## Progress Ledger

- Completed: strict native FastOMA conversion and full retained-row audit.
- Running at scheduler check: original QfO factorial reconciliation
  `21671_3`, assessment `21711`, corrected OrthoHMM `21706_0`, Proteinortho
  `21708`, SonicParanoid `21710`, and dedicated DGX scaling `21656_14`.
- Pending: corrected OrthoFinder `21706_1`, corrected legacy BLAST `21713`
  (resources), dependent factorial admission/conversion, and later scaling
  replicates. No unrelated jobs were stopped.
- Unresolved: original TreeFam-A trees/mapping remain unavailable; no
  family-level TreeFam uncertainty is claimed.
- Next: admit completed native runs before scoring; derive FastOMA's
  corrected-release tree from admitted corrected OrthoFinder output; run
  the frozen eight-cell uncertainty analysis only after all cells pass.

The overall publication goal remains incomplete.
