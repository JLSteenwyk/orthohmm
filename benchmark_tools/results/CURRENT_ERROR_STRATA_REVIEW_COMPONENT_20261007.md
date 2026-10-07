# Current Error-Stratum Main-Text Review Component

This new, local component exports the reviewed 7 October v3 manuscript and
direct assets from commit `3aaf36a4ed47384a48d50c624de21a9608a17f46`.
It does not replace rc5 or establish a complete executable study release.

## Actual Evidence

- [Fresh manuscript](PUBLICATION_MAIN_TEXT_20261007_v3.md),
  [21-page PDF](publication_main_print_20261007_v3/document.pdf) and
  [separate actual visual review](publication_main_pdf_visual_review_20261007_v3.json).
- [Build receipt](current_error_strata_review_build_20261007_v1.json):
  149 payloads, 9,272,442 payload bytes, 114 direct targets and 117 local links.
- [New archive](current_error_strata_direct_review_20261007_v1.tar.gz):
  6,354,070 bytes, SHA-256
  `42c8d5daced8cef058e1f31cd23f80476522be75184b3ef23ea17951d81b3a48`.
- `REVIEW_INDEX.json`: 54,171 bytes, SHA-256
  `41d31073e9c563ef0c530187d666a0b7e83d63d894038f19d54f2614877b8d53`.
- [Actual restoration](current_error_strata_review_restore_20261007_v1.json)
  checked all 150 regular members before extracting outside the checkout.
- [Copied verifier execution](current_error_strata_review_verify_20261007_v1.json)
  checked all 149 payloads under Python 3.10.13 with `-I -S -B`; exit 0,
  0.10 seconds, 24,576 KiB peak RSS and zero swaps in the
  [execution measurement](current_error_strata_review_verify_20261007_v1.time.txt).
  These are postprocessing costs, not native tool timings.

The new main text removes only the live progress-ledger hyperlink and adds
a version note. The [v2 failure](publication_main_pdf_review_failure_20261007_v2.json)
remains failed: changing the operational ledger invalidated its render-time
asset binding. No old receipt, PDF, source or archive was rewritten.
The [exact revision receipt](publication_main_review_asset_revision_20261007_v3.json)
and focused tests preserve every numerical table and scientific statement.

## Reproduction

Use a fresh destination and the two external digests above. The existing
standard-library [restorer](../restore_direct_review_archive.py) validates the
archive/index and every regular payload before extraction. Then actually run
the copied `benchmark_tools/bundle_publication_review.py`:

```sh
python3 -I -S -B /fresh/destination/benchmark_tools/bundle_publication_review.py \
  verify /fresh/destination \
  --manifest-sha256 41d31073e9c563ef0c530187d666a0b7e83d63d894038f19d54f2614877b8d53
```

HTML entrypoint: `benchmark_tools/results/publication_main_review_20261007_v3.html`.
Relative direct links are the portable navigation interface. PDF annotations
may retain historical workstation paths; those paths are provenance, not
portable link targets.

## Scope

The component contains the current fixed-bin figure, model-distance result
tables and its small retained evidence archive. It does not include every
linked document's transitive inputs, full raw predictions, native runtime,
or independent reexecution of feature inference, scoring or bootstrap draws.
Its automatic PDF-bounds receipt and all 21 raster images are included; the
separate manual visual receipt is retained in Git, not silently added to the
direct-link index. No OS containment, cross-host inference, biological
generalization, redistribution clearance, archival DOI or publication readiness
is established. Journal-specific presentation and scientific gaps remain open.
