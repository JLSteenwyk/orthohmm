# Corrected VGNC Reference-Block Mapping

Mapped all eight corrected-QfO native scored tables onto one reference-defined
block inventory. Every method reproduces its admitted precision, recall and
F1 to tolerance 1e-12, with exactly 23,934 TP+FN reference pairs. The
[generated count table](CORRECTED_VGNC_BLOCK_TABLE_20260926.md) preserves the
eight methods and their prediction semantics; no comparator was omitted.

The [machine-readable report](corrected_vgnc_blocks_20260926.json) has SHA-256
`01a51a62a536bcf759f2a0ef83c56196a92022b0bb8290448eacbfefa47d6005`.
Thirty-nine report/input/output record entries were rechecked. An independent
TSV readback sums each sparse table and reproduces all eight count triples.
Twenty-six related tests pass.

## Mapping and Checks

The frozen VGNC reference has 16,863 labels and 36,986 proteins. Eleven
proteins have multiple labels; merging labels that share reference proteins
produces 16,844 blocks. No prediction edges or scores determine those blocks.
All eight read-only prediction databases provide identical selected reference
protein/accession/species mappings, including the native last-label/alias
assignment. Each raw row is checked against that mapping and native eligibility
rules. TP and FN partition the mapped reference exactly. There are no observed
TP/FP or FN/FP overlaps in these eight tables.

The sparse tables aggregate TP/FP/FN by unordered block pair. Within-block
counts occupy diagonal cells. Cross-block rows are all FP; each contributes
once, not once per endpoint. Omitted cells have zero count contributions on
the complete fixed block grid, but are not claimed to have equal biological
eligibility. The reference inventory records protein count, species and
asserted-pair count per block. Both the reference inventory and eight sparse
tables remain local in `benchmarks/work/corrected_vgnc_blocks_20260926/`.

## Scientific Implications

Phylogenetic OrthoHMM has 19,660 TP, 13 FP and 4,274 FN; full OrthoFinder has
23,519 TP, 130 FP and 415 FN. This preserves the observed high-precision,
lower-recall trade-off rather than inferring superiority from FP alone. All
13 and 130 FP, respectively, cross reference blocks. High-sensitivity
OrthoHMM has 16,076 cross-block FP plus two within-block FP; the sequence-only
OrthoFinder checkpoint has 286,002 cross-block FP plus two within-block FP.
Within-block FP are retained: overlap-based block merging does not change the
native scorer's family labels or redefine its truth criterion.

Of the 16,844 blocks, 14,075 contain two proteins and one asserted pair.
Protein counts range from two to eleven, and asserted-pair counts from one
to fifteen. Only eight nonzero cross-block FP cells occur in phylogenetic
OrthoHMM, versus 38 in full OrthoFinder. This is substantially different from
the earlier [synthetic screen](DYADIC_F1_VARIANCE_RESULT_20260926.md), which
assigned 20 true pairs to every family. Its coverage cannot be transferred
to this sparse, imbalanced setting without further validation. Neither
protein disjointness nor these count tables establishes biological independence.

## Reproduction and Boundaries

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m benchmark_tools.map_corrected_vgnc_blocks --repo . \
  --output benchmarks/work/corrected_vgnc_blocks_20260926
python -m pytest -q tests/unit/test_map_corrected_vgnc_blocks.py \
  tests/unit/test_audit_qfo_vgnc_mapping.py \
  tests/unit/test_audit_vgnc_family_dependencies.py
```

Use a fresh destination; existing paths are refused. This is an aggregate-
bound raw-table readback, not a new native scorer execution or a requery of
all predicted edges. The pinned comparator manifest binds admissions and
aggregate JSON hashes. Supplemental raw hashes are established now, not
retrospectively asserted to have been checked by the historical admission.
Database checks cover selected reference mapping digests and size/mtime
stability during that query, not full database hashes. These limitations are
explicit in the report and must not be promoted to historical input proof.

No confidence intervals are calculated. A next uncertainty protocol must
address unequal reference sizes, rare false positives, near-boundary ratios,
native eligibility and possible dependence beyond shared families. Simulation
settings based on this fixed reference inventory can be label-independent;
choosing a correction to favor a method's observed errors cannot. Existing
point estimates, scientific defaults and the publication-readiness status are
unchanged. DGX remains deferred; no inference or timing run was launched.
