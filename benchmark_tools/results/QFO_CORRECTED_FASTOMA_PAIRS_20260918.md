# Corrected FastOMA Pair Preparation

Job **21742** is queued `afterok:21741`, the independent native-admission job.
It reserves 2 CPUs and 64 GiB for four hours on bizon, without requeue.
Actual corrected conversion has not yet run.

The frozen executor is
`benchmarks/work/publication_qfo_corrected_fastoma_pairs_v1` at
`6616e3a7ec4c46962ded0d34ce9b4150073a2210`. The retained batch script is
`qfo_corrected_fastoma_pairs_batch_20260918.sh`.

Preparation requires successful native-admission accounting and the exact
frozen native-admission source. It rechecks native inputs, outputs, helpers
and the corrected 78-species/984,137-accession universe before conversion.
It uses the same frozen QfO mapping as the other corrected comparators.

The converter consumes only `orthologs.tsv.gz`, never HOG clique expansion.
Every row must contain known cross-species input accessions. A temporary
disk-backed index canonicalizes orientation, removes duplicates and records
native row, distinct pair and duplicate relation counts separately. The native
row count must equal the independent admission count. The entire gzip is read
before any pairs are emitted by the distinct converter.

Mapping loss or count mismatch fails preparation and preserves partial evidence.
Checksummed sources are verified again before successful final files are
published. There is no implicit restart, overwrite or silent score transfer.

Prospective outputs under
`benchmarks/results/qfo_corrected_comparator_pairs_v1/fastoma/`:

- `pairs.tsv`: sorted distinct native pairs.
- `pairs.qfo.tsv`: reference-filtered pairs; successful corrected preparation
  requires no loss and therefore identical pair content.
- `results.json`: counts, provenance, mapping and admission accounting; status
  `corrected_fastoma_pairs_prepared_unscored`, accuracy and publication flags false.

The manifest explicitly records that these are native phylogenetically inferred
pairs using a supplied corrected OrthoFinder species tree. It does not imply
independent FastOMA tree inference or comparable dedicated runtime.

All 45 focused pair-preparation/deduplication tests and Bash syntax validation
pass. Tests cover admission scope/source-record checks, duplicate handling,
count drift, malformed/foreign rows, mapping loss and no overwrite. The earlier
complete historical distinct-pair audit provides real-data evidence for the
underlying converter, not an actual corrected conversion result.

Next: extend the frozen six-endpoint scorer and independent score admission for
this manifest type, then queue them after successful conversion. No endpoint,
mapping or historical score has changed.
