# Revised Main-Text Review

The [current six-page PDF](publication_main_print_20260930_v1/document.pdf)
and [citation-rendered HTML](publication_main_review_20260930_v1.html)
include the restored-assets full OrthoBench reproduction and GO/EC composition
analysis added to the main Markdown on 30 September. All six PDF pages were
visually inspected; the bounds check found no violations. The
[visual-review receipt](publication_main_visual_review_20260930_v1.json)
records the scope: main-text pages, not a repeated review of every linked figure
or validation of scientific claims. Figures remain linked, not embedded.

The [actual archive/fresh-extraction receipt](publication_review_component_20260930.json)
verifies 45 payload files, 3,722,238 bytes, all 25 direct local targets and
26 HTML link occurrences. Review artifacts and workflow use committed revision
`4b7c35519a19d3f23ec26d8570899b7203e5392e`. The main text and render-time ledger
are bound to `52b93750b0596a1702763a6b6fb356b1953bc571`; the later ledger is
not substituted into historical rendering evidence.

Local archive:
`benchmarks/work/publication_review_component_20260930/orthohmm-main-review-20260930-4b7c3551.tar.gz`

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Archive | 2,134,463 | `e08112c571cec6ab54fea6a896fa11a8f3e32b020fe35d13b9f17929eb8ab902` |
| `REVIEW_INDEX.json` | 17,136 | `64a9530090d4b016ef9d692d2ad96b80734a1c21885e6f8e78b7bdc245220ffd` |

The self-generated archive was checked for regular/directory-only safe member
paths, extracted to a fresh temporary directory and verified using
`/usr/bin/python3 -I -B` with PATH `/no-git`. Its complete result equals the
original component verification. Temporary extraction was removed; the archive
and original component remain retained. No checkout, Git, Pandoc, browser or
scientific Python packages are needed for this verification.

## Reproduce

The [component guide](../PUBLICATION_REVIEW_COMPONENT.md) now accepts explicit
render/print/review receipt paths. The selected stages are stored in a schema-v2
index and their render/PDF binding is checked. The historical default remains
schema v1. 74 focused exporter/render/print/PDF-review tests pass, including
explicit selection, malformed mappings, missing chain bindings and isolated
verification after removing a synthetic source repository.

The [actual historical-component check](publication_review_legacy_compatibility_20260930.json)
also proves that the updated verifier accepts the retained 43-file v1 component,
with a result exactly equal to its original admission. Its inference, rendering
and archive extraction were not repeated.

After transferring and extracting the new component:

```sh
python3 -I -B /relocated/component/benchmark_tools/bundle_publication_review.py \
  verify /relocated/component \
  --manifest-sha256 64a9530090d4b016ef9d692d2ad96b80734a1c21885e6f8e78b7bdc245220ffd
```

Use the relative HTML entrypoint for portable direct-link navigation. The
preserved PDF may retain workstation-specific link annotations. The new review
does not replace the older archive or modify its receipts.

## Boundary

This is a local main-text review component, not full transitive study evidence,
an executable all-method release, cross-host inference, rights clearance or
public archival deposition. Native scores, methods and timing identities are
unchanged. The main draft still lacks controlled timing and appropriate
comparison uncertainty for several QfO endpoints; final journal formatting
and submission readiness remain open.

A [fresh read-only host observation](threadripper_preflight_20260930.json)
found about 80.99 competing CPU-core equivalents despite an empty Slurm queue,
including unrelated IQ-TREE and Python analyses. One process-read error was
retained. This is a three-second, non-atomic observation, not a whole-run
certificate or a forecast. No process was signalled and no production timing
or overhead experiment was launched. Quiet-window coordination was requested.
