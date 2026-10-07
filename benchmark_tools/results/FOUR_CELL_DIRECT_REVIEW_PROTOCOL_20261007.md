# Four Cell Main Review Export Protocol

Export the already reviewed v4 manuscript and its direct assets using the
unchanged existing exporter. The older v3 review archive does not contain this
new four-cell source. This is remaining manuscript-delivery work, not new
science, repetition of an old restore or a full-study release claim.

## Selection

Review, workflow and ledger-selection revision:
`3b1dc7add7e8f118fbb653c0b65866c17a65111d`.
There is no live-ledger hyperlink in this source; selecting that revision does
not add the changing operational ledger to the direct-asset inventory.

- Main: `benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261007_v4.md`,
  SHA256 `0b7012e3dd46bd3b028857bdb3cbe3a9e9048936b90f3df9e4833cb342d5705a`.
- Render: `benchmark_tools/results/publication_main_review_20261007_v4_checked_assets.json`,
  SHA256 `995f707cf3dcf43b73bc6f95a8fbdfb12decb10884e0ed11142a650e7927fb21`.
- Print: `benchmark_tools/results/publication_main_review_20261007_v4_checked_print/print.json`,
  SHA256 `e4abfea87036c3e4252267603bcf7166d766ebfbb13f5c712e23d4dd67f0d681`.
- PDF review: `benchmark_tools/results/publication_main_review_20261007_v4_checked_pdf_review/report.json`,
  SHA256 `88e054a69db01f6ae73e06ed78c641d6036d5d78a0bc7a136990d37ed9de5a8d`.
- Exporter: `benchmark_tools/bundle_publication_review.py`, unchanged SHA256
  `4b9646bf7fd88c750bf9df4e7fad22c4fac6d5f9810cf211a639c15009a07999`.
- Restorer: `benchmark_tools/restore_direct_review_archive.py`, unchanged SHA256
  `e14d4ee8fd9f2cafbda64ca2d2f2ae964611df16734c6bc113de3070b1e89090`.

## Execution

Commit/push this selection before one selected export. Use the existing clean
environment and scientific Python 3.10 with `-I -S -B`:

```sh
python -I -S -B benchmark_tools/bundle_publication_review.py build \
  --repo . \
  --review-revision 3b1dc7add7e8f118fbb653c0b65866c17a65111d \
  --ledger-revision 3b1dc7add7e8f118fbb653c0b65866c17a65111d \
  --workflow-revision 3b1dc7add7e8f118fbb653c0b65866c17a65111d \
  --main-text benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261007_v4.md \
  --render-receipt benchmark_tools/results/publication_main_review_20261007_v4_checked_assets.json \
  --print-receipt benchmark_tools/results/publication_main_review_20261007_v4_checked_print/print.json \
  --review-receipt benchmark_tools/results/publication_main_review_20261007_v4_checked_pdf_review/report.json \
  --output benchmarks/work/four_cell_direct_review_20261007_v1
```

Build from committed Git blobs only. Expect the exact v4 entrypoint and 107
direct targets/110 local links/21 PDF pages; the checked index determines the
complete regular-payload count and digest. No guessed digest is admissible.

After success, feed only the checked index's regular file names plus
`REVIEW_INDEX.json` to GNU tar using NUL-delimited `--files-from=-`,
`--no-recursion`, fixed mtime/owner/group and a fresh destination
`benchmark_tools/results/four_cell_direct_review_20261007_v1.tar.gz`.
Use structured JSON parsing, not a recursive directory tar that would include
unexpected directories or browser state. Retain the actual archive digest and
command, then commit/push the successful build/archive receipts and index before
using those external anchors for a fresh restore.

Existing restorer must check the entire archive before extraction into fresh
`/tmp/orthohmm_four_cell_review_20261007_v1`, outside the checkout. Retain its
receipt at `benchmark_tools/results/four_cell_direct_review_restore_20261007_v1.json`.
Byte-check the copied verifier and all payloads, then actually invoke that
copied file with `-I -S -B`, the externally retained index digest and `/tmp` as
working directory. Record actual exit/stdout/time; checks that precede copied
execution do not prove copied execution. Preserve any failure without an
automatic retry or editing its sources, results or inventory.

## Boundaries

The exporter includes exact direct HTML targets, the source, bibliography,
PDF, all 21 page rasters, checked stage receipts and existing standard-library
support. The separate later manual visual record and new claim addendum remain
in Git, outside this unchanged direct-link export contract. Never silently add
extra files to its index or claim that it includes transitive inputs.
The copied verifier checks preserved bytes and direct links; it does not rerun
scoring, plots, rendering, inference or biological validation. Whole-study
raw-data/runtime closure, public redistribution, journal formatting, external
deposition and scientific limitations remain distinct. Older v3/rc5 bundles and
all ongoing native jobs remain untouched. No archival DOI, isolated speedup or
publication readiness follows from this export.
