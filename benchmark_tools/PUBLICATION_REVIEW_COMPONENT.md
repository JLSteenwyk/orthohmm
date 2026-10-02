# Main-Text Review Component

`bundle_publication_review.py` exports the retained main-text Markdown, HTML,
six-page PDF, bibliography, direct linked documents/figures and render/print/PDF
receipts with page images. It reads committed Git blobs, not current worktree
files, and preserves their relative paths. The independently pinned index
permits offline verification using only the Python standard library.

For the later seven-page review, use the explicit committed stage selection
and index digest in the [2 October archive receipt](results/PUBLICATION_REVIEW_COMPONENT_RESULT_20261002.md).
That 55-file snapshot passes fresh-extraction verification outside the checkout;
the historical defaults and older archives below remain unchanged.

The default historical review revision is `e2e96746b21d552ea46b9ce63634519f9bae1f82`. The ledger
snapshot used during rendering is its parent,
`fdfe7fb70c583a6ad4ea5c4ae06dd2663d31b07f`; the later ledger must not replace
it just because it has more current status. Both revisions are checked against
the retained render/print/PDF identities. Select a committed workflow revision
containing this exporter, this guide and the project license:

```sh
python -B benchmark_tools/bundle_publication_review.py build \
  --repo . --workflow-revision COMMITTED_WORKFLOW_REVISION \
  --output /absolute/fresh/main-text-review
```

For a newer committed review, explicitly select all three repository-relative
stage receipts. Do not rename or overwrite old receipts. The ledger revision
must supply the exact snapshot pinned during rendering, which may precede
the commit that added the new review artifacts:

```sh
python -B benchmark_tools/bundle_publication_review.py build \
  --repo . --review-revision NEW_REVIEW_COMMIT \
  --ledger-revision ACTUAL_RENDER_TIME_LEDGER_COMMIT \
  --workflow-revision COMMITTED_WORKFLOW_REVISION \
  --render-receipt benchmark_tools/results/publication_main_render_20260930_v1.json \
  --print-receipt benchmark_tools/results/publication_main_print_20260930_v1/print.json \
  --review-receipt benchmark_tools/results/publication_main_pdf_review_20260930_v1/report.json \
  --output /absolute/fresh/new-main-text-review
```

Explicit selection produces a `publication_direct_review_v2` index containing
the stage-path mapping and verifies the selected render/print/PDF chain.
The default historical export remains schema v1; both formats are supported
by the current offline verifier. Neither selects files from a dirty worktree.

The result reports the index SHA-256 and relative HTML/Markdown/PDF entrypoints.
After relocation, supply that external digest, rather than trusting a digest
read from the bundle itself:

```sh
python3 -I -B /relocated/main-text-review/benchmark_tools/bundle_publication_review.py \
  verify /relocated/main-text-review --manifest-sha256 EXTERNAL_INDEX_SHA256
```

Open the relative HTML entrypoint to follow its direct local links. The PDF is
preserved byte-for-byte and may retain original workstation-specific link
annotations; its links are not the portable navigation interface. Source
receipts also retain original absolute paths as provenance, not active reads.
The verifier needs no Git, checkout, Pandoc, Chrome, native tools or packages.

## Boundary

This is a main-text review component, not a complete study archive or a public
release. Local links inside linked Markdown documents, transitive prediction
tables and raw datasets are not included. External URLs, anchors, scientific
contents and redistribution rights are not certified. Project license inclusion
does not clear third-party assets. Exporting does not rerun native inference,
scoring, plotting, rendering or visual review, and cannot establish independent
generalization, controlled timing or publication readiness. Historical negative
results and explicit limitations remain unchanged.
