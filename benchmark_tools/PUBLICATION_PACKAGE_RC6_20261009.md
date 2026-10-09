# OrthoHMM Study Candidate 2026.10.09-rc6

This local candidate preserves every indexed rc5 payload from its immutable
package. Its former root README, selection, verifier and index are retained
under `history/rc5/`. The scientific revision is unchanged. All inherited
review routes and their embedded status statements are historical; this guide
selects the new reporting revision rather than editing those snapshots.

## Current Reporting

The current [22-page manuscript](terminal-direct-review/benchmark_tools/results/native_qfo_terminal_failures_20261009_v1_print/document.pdf),
[HTML review](terminal-direct-review/benchmark_tools/results/native_qfo_terminal_failures_20261009_v1_review.html)
and [Markdown source](terminal-direct-review/benchmark_tools/results/native_qfo_terminal_failures_20261009_v1_manuscript.md)
retain four admitted native QfO cells and three missing cells. Native11 failed
scoring with OUT_OF_MEMORY; native12 inference exited -11 (SIGSEGV). No failed
result was scored or retried. The four admitted scores, uncertainty intervals
and retained figures are unchanged. A complete fresh factorial remains unavailable.

`terminal-direct-review/` contains exactly 146 payloads plus its anchored index,
including 110 direct local targets and 114 local HTML link occurrences. Its
HTML links resolve locally. Links inside linked documents, historical PDF
annotations and absolute provenance paths are not guaranteed portable. The
terminal-review manifest contains file identities and outcome/resource metadata,
not raw sequences or expanded process streams.

Separate `evidence/rc6/results/` records contain the actual manual page review,
content readback, citation inventory, TSV exports and archive/restoration/copied
verification receipts. They are not silently added to the strict direct index.
Those dated repository-evidence documents retain their original relative links
and pre-execution status text; use the current entrypoints above and subsequent
execution receipts. Workflows/tests under `evidence/rc6/` preserve repository-layout
provenance, not a standalone complete scientific environment.

## Restore And Verify

Keep archive and index digests independently. With a trusted, byte-checked
standard-library package reader and a fresh destination:

```sh
python3 -I -S -B /trusted/bundle_publication_package.py restore /candidate.tar.gz \
  --archive-sha256 ARCHIVE_SHA256 --manifest-sha256 PACKAGE_INDEX_SHA256 \
  --output /fresh/package
python3 -I -S -B /fresh/package/bundle_publication_package.py verify /fresh/package \
  --manifest-sha256 PACKAGE_INDEX_SHA256
```

Outer verification checks selected bytes and inventories. It does not execute
nested components or reproduce scientific analyses. The separate direct
component can be checked with its copied, byte-checked verifier:

```sh
python3 -I -S -B /fresh/package/terminal-direct-review/benchmark_tools/bundle_publication_review.py \
  verify /fresh/package/terminal-direct-review \
  --manifest-sha256 506c2a0b64d95282debb4eeae5b111fcaae39dbf8f5ebeb5e96f0319a227e396
```

The pre-package direct restoration and actual copied-verifier execution are
retained receipts, not proof that this outer candidate has executed that check.
Consult its dated execution records for what actually ran. No instruction here
authorizes inference retries, new scoring/admission or repeated source searches.

## Scientific Limits

QfO and OrthoBench are development-exposed primary benchmarks; Three Kingdoms
is supplementary. The six-endpoint mean is a project-defined secondary summary,
not official QfO F1. SwissTrees adjusted intervals include zero; partial native
factorials do not establish a general HMM or phylogeny advantage. Independent
clade evidence is not proof of family-disjoint confirmation. The original
TreeFam-A reference reconstruction and several uncertainty endpoints remain
unresolved. No new default or superiority claim follows from this package.

Timing measurements were collected on a shared Threadripper while other
analyses were running. Competition for CPU, memory bandwidth and I/O may have
affected elapsed times, with an unknown and potentially tool-dependent impact.
These are observed shared-host timings, not estimates of isolated performance.
No DGX, quiet host or interference with unrelated jobs is required.

Complete transitive raw-data/runtime restoration, all-method native reproduction,
redistribution clearance and public deposition remain unproved. This is a local
versioned reporting candidate, not a hermetic full-study archive, publication-ready
release, archival DOI or journal submission. Completion requires the full
seven-part scientific requirement audit, not successful byte verification alone.
