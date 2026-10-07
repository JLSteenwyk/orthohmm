# Four Cell Manuscript Review Archive And Restoration

The exact reviewed v4 manuscript and direct assets now have an independently
anchored archive that restored and verified outside the checkout. This closes
the current manuscript's direct-review delivery gap; it does not complete the
whole-study executable package or establish publication readiness.

## Actual Evidence

| Step | Retained Evidence | Actual Result |
| --- | --- | --- |
| Export from committed review | [Build receipt](four_cell_direct_review_build_20261007_v1.json), [prospective selection](FOUR_CELL_DIRECT_REVIEW_PROTOCOL_20261007.md) | Unchanged exporter, exact review commit `3b1dc7add7e8f118fbb653c0b65866c17a65111d`; 142 payloads, 8,916,531 bytes, 107 direct targets, 110 local links and 21 PDF pages; exit 0 |
| Regular-file archive | [Archive](four_cell_direct_review_20261007_v1.tar.gz), [actual command receipt](four_cell_direct_review_archive_20261007_v1.json), [external index copy](four_cell_direct_review_index_20261007_v1.json) | 143 logical regular members including the index; archive 6,158,614 bytes; no recursive browser-state inclusion; exit 0 |
| Fresh external restoration | [Restorer receipt](four_cell_direct_review_restore_20261007_v1.json) | All archive/index/payload digests and modes checked before extraction to `/tmp/orthohmm_four_cell_review_20261007_v1`; exit 0 |
| Actual copied verification | [Execution receipt](four_cell_direct_review_verify_20261007_v1.json) | Copied verifier byte-checked, then executed with `-I -S -B` and `/tmp` as working directory; all 142 payloads and direct links match; exit 0 |

Protocol `ccec7e4b` was committed/pushed before export. Actual archive, index
and build/archive receipts were committed/pushed as `0a890b7e` before external
restore. No historical source, component, archive, manuscript, PDF or scientific
result was modified. Separate later manual visual review and the current
four-cell claim addendum remain in Git, not silently added to this index.

External anchors:

- Archive SHA256: `ac604c4d16ec4bb7aead2f75adb43ad8364e6b0d9d359071b6f66a94bea8c59e`.
- `REVIEW_INDEX.json` SHA256: `3a6e46bc6aafea5ba3b012e78091b778d046fe43790fa9cf0c1413eeeec9a4a1`,
  52,225 bytes.
- Copied verifier SHA256: `4b9646bf7fd88c750bf9df4e7fad22c4fac6d5f9810cf211a639c15009a07999`.

Measured postprocessing costs, not native tool timings: export 1.87 seconds
and 23,040 KiB maximum RSS; archive wrapper 0.61 seconds and 23,828 KiB;
restore 0.22 seconds and 15,360 KiB; copied verification 0.10 seconds and
23,040 KiB. Actual copied result equals the original export result, including
all entrypoints and scope flags.

## Reproduction

Use the unchanged [standard-library restorer](../restore_direct_review_archive.py)
and a fresh destination. Retain the digests above externally rather than
trusting an index supplied inside an unverified archive:

```sh
python3 -I -S -B benchmark_tools/restore_direct_review_archive.py \
  benchmark_tools/results/four_cell_direct_review_20261007_v1.tar.gz \
  --archive-sha256 ac604c4d16ec4bb7aead2f75adb43ad8364e6b0d9d359071b6f66a94bea8c59e \
  --index-sha256 3a6e46bc6aafea5ba3b012e78091b778d046fe43790fa9cf0c1413eeeec9a4a1 \
  --output /fresh/destination --receipt /fresh/restore-receipt.json

python3 -I -S -B /fresh/destination/benchmark_tools/bundle_publication_review.py \
  verify /fresh/destination \
  --manifest-sha256 3a6e46bc6aafea5ba3b012e78091b778d046fe43790fa9cf0c1413eeeec9a4a1
```

HTML entrypoint: `benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261007_v4_checked.html`.
Portable navigation uses relative HTML links. The preserved PDF's annotations
may retain workstation paths; those are historical provenance, not portable
links. The index and verifier retain original absolute source paths as metadata
while checking copied payloads at the new location.

## Remaining Scope

This is preserved reporting-content verification, not raw prediction admission,
feature inference, scoring, bootstrap or plot/render reexecution. Linked
documents' transitive inputs, complete native tool/runtime closure and some
raw provenance remain separate. No original-checkout file-access trace, OS
containment, cross-host execution or third-party redistribution clearance is
established. The existing method's scientific limitations, pending native
outcomes, other-endpoint uncertainty and archival/submission steps remain open.
No DOI, isolated performance ranking, default promotion or publication
readiness is claimed. The full publication goal remains active and incomplete.
