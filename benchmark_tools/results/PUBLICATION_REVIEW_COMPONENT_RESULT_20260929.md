# Relocated Main-Text Review Component

The [executed receipt](publication_review_component_20260929.json) records an
actual local export, archive and fresh-extraction verification. This closes
the main review's direct-file checkout dependency, not the full study archive.

| Artifact | Verified retained scope |
| --- | --- |
| Main text | Exact Markdown and citation-rendered HTML |
| Printed review | Exact six-page PDF and six page images with print/bounds receipts |
| Direct local targets | All 23 unique targets, 24 HTML occurrences, including three figure PDFs, comparison table and linked result/protocol documents |
| Support | Bibliography, render/print/review source, exporter, guide and project license |

There are 43 payload files totaling 3,581,452 bytes. Review artifacts use
`e2e96746b21d552ea46b9ce63634519f9bae1f82`; the ledger snapshot that was
actually read during rendering uses parent
`fdfe7fb70c583a6ad4ea5c4ae06dd2663d31b07f`. Exporter/guide/license use
`73ddcdbb` (full revision in the receipt). The newer ledger is deliberately
not substituted into historical render evidence.

Local archive:
`benchmarks/work/publication_review_component_20260929/orthohmm-main-review-e2e96746-workflow-73ddcdbb.tar.gz`

- Archive bytes: 2,031,512.
- Archive SHA-256: `9a29bfb4cc9095c3063bb9dd7ab6e5790b1824225671327c872e3f7fa4ed318b`.
- `REVIEW_INDEX.json` bytes: 16,130.
- Index SHA-256: `d6cc0e930c61ec9b97a8124284a9974ec9176e626474f78941b18dc1a346c7f5`.

The archive was extracted to a fresh temporary directory and the bundled
verifier ran under `/usr/bin/python3 -I -B` with PATH `/no-git`. Its full
verification result equals the original result. Temporary extraction was
removed afterward; the original component and archive remain retained.
Historical absolute paths are not read, and Pandoc/Chrome/project dependencies
are not needed for verification.

## Reproduce

See the [component guide](../PUBLICATION_REVIEW_COMPONENT.md) for the build
command and exact historical revisions. To verify after transfer/extraction:

```sh
python3 -I -B /relocated/component/benchmark_tools/bundle_publication_review.py \
  verify /relocated/component \
  --manifest-sha256 d6cc0e930c61ec9b97a8124284a9974ec9176e626474f78941b18dc1a346c7f5
```

Open `benchmark_tools/results/publication_main_review_20260929_v3.html`
inside the extracted component for portable direct-link navigation. The PDF
is preserved byte-for-byte and may retain workstation-specific link annotations.

61 focused exporter/render/print/PDF-review tests pass. New component tests
include deletion of the source repository before isolated verification, dirty
checkout exclusion, wrong historical ledger selection, repeat refusal, changed,
missing, extra and symlink files, mode differences, external-index anchoring,
inconsistent stage evidence, escaping URLs and rendered-link inventory mismatch.
Synthetic test PDFs/images are integrity fixtures, not visual-validation evidence.

## Boundary

No native inference, scoring, plotting or visual review was repeated. The
component preserves the existing draft and limitations. Links inside linked
documents, raw data, prediction tables and full transitive evidence are not
included. External URLs, fragment anchors and scientific claims are not audited
by file hashing. Project license inclusion does not clear third-party assets.
This is a local review archive, not public deposition, a versioned package
release, cross-host restoration, controlled timing or publication readiness.
All 823 source pins for deferred audit 22379 still match; none was edited.
