# Main-Text Review Component

`bundle_publication_review.py` exports the retained main-text Markdown, HTML,
six-page PDF, bibliography, direct linked documents/figures and render/print/PDF
receipts with page images. It reads committed Git blobs, not current worktree
files, and preserves their relative paths. The independently pinned index
permits offline verification using only the Python standard library.

The review revision is `e2e96746b21d552ea46b9ce63634519f9bae1f82`. The ledger
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
