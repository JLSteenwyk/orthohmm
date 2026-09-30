# Current Main-Text Review Snapshot

The [six-page PDF](publication_main_print_20260930_v2/document.pdf) and
[citation-rendered HTML](publication_main_review_20260930_v2.html) now include
the exact frozen scoring specification, completed eight-method FAS coverage
audit and all-listed-tag-name TreeFam context search. They preserve native
benchmark scores, negative findings and unresolved timing/uncertainty claims.
All six pages were visually inspected; the bounds check found zero violations.
The [manual review receipt](publication_main_visual_review_20260930_v2.json)
records readable body text, headings and bibliography without observed clipping
or overlap. Figures remain linked, not embedded or newly visually revalidated.

An [initial selector check](publication_main_pdf_selector_failure_20260930_v2.json)
failed because its requested phrase "scoring specification" differs from the
literal source wording "numerical specification". The failed check is retained.
The separately retained corrected check uses the literal wording without
changing or reprinting the manuscript/HTML/PDF; all four selectors appear.
This is an operator selector correction, not a changed scientific criterion.

## Archive And Reproduction

The [archive/fresh-extraction receipt](publication_review_component_20260930_v2.json)
verifies **48 payload files**, 3,892,072 bytes, 28 direct local targets and
29 HTML link occurrences. Artifacts/workflow use committed revision
`e23faf23e63f40a0e060cf29cb1cca546da39234`. Main text and render-time ledger
are bound to `706786b12039ac1c66ef23d3104f88190ec86a51`; the newer ledger is
not substituted into the historical render inputs.

Local archive:
`benchmarks/work/publication_review_component_20260930_v2/orthohmm-main-review-20260930-e23faf23.tar.gz`

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Archive | 2,257,011 | `c84ba6bffcf317ef6912c15d743fd68b06a247f26144e0a6663f524c570a927c` |
| `REVIEW_INDEX.json` | 18,187 | `fda371438141ffc0d85c271e48457a3439375c904ac919451efaaa5df2b13c5e` |
| PDF | 124,470 | `427ddd55bfa91e039565a1841d5d7163fa319b69cb7c22547ab85fe3757e610d` |

All archive members were checked for safe, unique regular-file/directory paths
before data-filtered extraction into a fresh temporary directory. Isolated
`/usr/bin/python3 -I -B` verification with PATH `/no-git` yields exactly the
same result after extraction. Temporary extraction is removed; the archive
and original component remain retained. No checkout, Git, Pandoc, browser or
scientific Python package is needed for verification. Older archives/receipts
are neither overwritten nor reinterpreted as current scientific results.

After transfer and extraction, verify with the external index digest:

```sh
python3 -I -B /relocated/component/benchmark_tools/bundle_publication_review.py \
  verify /relocated/component \
  --manifest-sha256 fda371438141ffc0d85c271e48457a3439375c904ac919451efaaa5df2b13c5e
```

Use the relative HTML entrypoint for direct-link navigation. The PDF preserves
workstation-specific link annotations. The unchanged exporter supports the
explicit render/print/review stage paths documented in the
[component guide](../PUBLICATION_REVIEW_COMPONENT.md); prior focused workflow
tests were not repeated merely because the goal resumed. Actual new build and
fresh-extraction verification check this new input snapshot directly.

## Boundary

This is a main-text/direct-asset review component, not a complete executable
study release, cross-host inference reproduction, rights clearance or public
archival deposition. Linked Markdown documents' transitive evidence is not
included. Separate manual-review and failed-selector receipts remain in Git,
not in this direct-asset payload; PDF bounds evidence is included.

No inference, scoring, plotting, completed calibration or source search was
repeated. Scientific code/defaults, native scores and timing identities are
unchanged. Controlled Threadripper timing remains deferred without a new host
probe, coordination question, unrelated job/service action or DGX access.
Other-QfO uncertainty, original TreeFam sources, high-CPM admission and final
release work remain open. This review advances packaging but does not satisfy
the full publication goal or establish submission readiness.
