# Native-Evidence Review Component Restored

The unchanged committed exporter builds the actual18-page main-text review
from source/review revision `c2a0ea998c04122f70213116d47d396ae247a048`, preserving
the exact render-time ledger at `47b7f14ba0c2e8c5b7f070f5b11620dd356a1eb9`.
This new component does not modify the older rc4 study candidate or require
another whole-study archive build while native inference remains active.

The [public index](native_main_review_component_index_20261006.json) is an exact
copy of canonical `REVIEW_INDEX.json`:41,307 bytes, SHA256
`6dffd552fa40a7a3f0f02fae8f291482264d6789773f9c4df096498b2808db14`.
It selects113 regular payloads totaling8,143,531 bytes, including85 direct
targets,88 HTML-link occurrences,18-page PDF,14 retained page images and the
selected render/print/bounds receipts. Native OrthoBench/QfO point snapshots,
guarded family-interval bindings/readbacks, native QfO vector figure and
functional-pair composition result/readback are directly included.
Each exported payload retains its exact Git revision/blob, bytes, digest and
mode. This is not a transitive raw prediction/reference/runtime archive.

## Actual Archive And Copied Verification

The local archive contains114 regular file members:113 selected payloads and
the external-anchor index. It is4,620,877 bytes, SHA256
`95fed303f76a40f281338c519be7a58304a752b1316c9c0714336d8a9dbe6df3`:

```text
benchmarks/work/native_main_review_component_20261006_v1.tar.gz
```

It remains local, not committed or uploaded. The
[actual execution receipt](native_main_review_component_execution_20261006.json)
records build/verification commands and results, copied archive/index/reader
identities, external restoration location and independently checked113 restored
payloads. Standard-library `tarfile` creates the archive with exclusive output;
the copied archive is hash-checked before extraction. Validate exact member
names, regular types, byte identities and payload modes before using the
Python3.12 data extraction filter. The completed archive is built only once.

The [initial staging refusal](native_main_component_initial_restore_refusal_20261006.json)
occurs before extraction or copied verification. Its inline guard incorrectly
applies the payload-only0644/0755 rule to `REVIEW_INDEX.json`, whose actual
mode is0664. The unchanged component verifier externally anchors this regular
index by digest, without imposing the payload-only mode rule. Preserve the
complete archive, copied bytes and original empty restoration directory;
resume extraction with payload modes checked against each manifest row and
the observed index mode explicitly accepted. Do not change any payload, index,
scientific source, hash rule or validator, and do not rebuild the archive.

Restoration and copied verification run outside the checkout, under
`/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/tmp/build_tmp/orthohmm-native-main-component-20261006-0o8t_a9j/restored`.
The copied verifier exits zero with `PATH=/no-git`, isolated/site-disabled
Python and no inherited project/scientific-package paths. It needs no Git,
Pandoc, Chrome, scientific packages, raw references or original checkout.
The [retained syscall trace](native_main_review_component_verify_20261006.strace)
is byte-identical to the original431,153-byte trace. It shows accesses to the
copied index/all113 payload paths and no original workspace or goal-attachment
prefix. This observation covers one verification command, not OS containment,
all-command certification or a bundled Python runtime. System-library reads
remain allowed. Copy commands emit GNU `cp -n` portability warnings, not failed
copy/analysis outcomes; the public index/trace hashes are independently checked.

## Reproduce And Navigate

Rebuild with a fresh destination from the recorded commit, never dirty worktree
bytes. Use the actual explicit schema3 selection, not the component guide's
dated historical default6/7-page examples:

```bash
python3 -I -S -B benchmark_tools/bundle_publication_review.py build \
  --repo . --review-revision c2a0ea99 --ledger-revision 47b7f14b \
  --workflow-revision c2a0ea99 \
  --main-text benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261006.md \
  --render-receipt benchmark_tools/results/publication_main_render_20261006_v1.json \
  --print-receipt benchmark_tools/results/publication_main_print_20261006_v1/print.json \
  --review-receipt benchmark_tools/results/publication_main_pdf_review_20261006_v1/report.json \
  --output /absolute/fresh/native-main-review
```

After trusted restoration, use the separately retained external anchor:

```bash
python3 -I -S -B /restored/benchmark_tools/bundle_publication_review.py \
  verify /restored \
  --manifest-sha256 6dffd552fa40a7a3f0f02fae8f291482264d6789773f9c4df096498b2808db14
```

Portable navigation entrypoint:
`benchmark_tools/results/publication_main_review_20261006_v1.html`.
The Markdown and PDF retain their original relative paths. Absolute paths in
JSON provenance and PDF link annotations are historical metadata, not verifier
reads; use the relative HTML links. Links inside linked Markdown documents
may remain unavailable. External URLs, anchors and redistribution rights are
not certified. The separate manual visual-review receipt is not a direct render
target and is not silently added to this immutable component; its bounded
14-page inspection remains in the repository's actual review result.

## Checks And Full-Goal Status

All52 joined cases pass in5.52s with no failure/error/skip, including seven new
actual-component cases. [JUnit](native_main_component_tests_20261006_v1.xml).
Independently compare all113 public-index payload identities with their actual
Git blobs, all114 archive members against selected bytes/modes, recorded
source/ledger/stage selection, actual copied-verifier output, trace accesses
and preservation of the initial refusal. Set `ORTHOHMM_NATIVE_REVIEW_ARCHIVE`
to the retained local archive to enable its byte-level archive test; absent
that explicit private/local artifact, only that case skips. The other checks
do not depend on the temporary restoration directory. Earlier19 source/five
PDF cases and explicit-main exporter fixtures also pass. These are reporting
and packaging checks, not new scientific admissions or generalization tests.

Fresh scheduler confirms22444 RUNNING2:21:32/native2:20:40; original22445 stays
dependency-pending. No inference/scoring/recount/bootstrap or timing retry
occurs. Observe those original handles; do not duplicate them or advance
index9 before actual terminal review and fresh capacity/accounting gates.
Remaining native QfO cells, other-endpoint uncertainty, provenance/independent
validation limitations and complete versioned executable study integration
remain in scope. There is no public upload, archival DOI, rights clearance or
publication-readiness claim. Shared-host timings continue to disclose unknown,
potentially tool-dependent CPU/memory-bandwidth/I/O contention; no quiet host
or DGX is required. The full publication goal stays active and incomplete.
