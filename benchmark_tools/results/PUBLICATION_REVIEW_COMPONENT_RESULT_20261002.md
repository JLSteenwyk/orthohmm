# Seven Page Manuscript Review Archive

The [reviewed seven-page main text](PUBLICATION_MAIN_REVIEW_20261001.md) now
has a portable direct-asset archive. It includes the completed parameter panel,
FAS mixture diagnostic and numerical-restoration scope already present in that
review. No manuscript rendering, figure review, inference or scoring is repeated.
Older archives, negative results and unresolved claims stay unchanged.

## Archive Verification

The [machine-readable receipt](publication_review_component_20261002.json)
records the actual build, archive and fresh-extraction verification. The archive
contains **55 payload files**, 4,262,679 bytes, 34 direct local targets and
35 HTML link occurrences. All payload identities, modes and render/print/review
bindings pass the unchanged verifier. The separate manual visual receipt remains
in Git; its seven page images and stage identities match the archived payload.
This reuses the recorded visual review, not a new visual inspection of figures.

Local archive:
`benchmarks/work/publication_review_component_20261002/orthohmm-main-review-20261002-c98b11fb.tar.gz`

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Archive | 2,526,021 | `23bcb9aa8920e1df642692d7d919f183aa4098fe1c4e22194856fd61722c20d4` |
| `REVIEW_INDEX.json` | 20,617 | `d2b696bdeb9d70f13e1d5b7ff6e52808c0267f5f2f28b4e56049b654e2d5557a` |
| PDF | 134,445 | `2330f990029a71542d1a47446d82c566264f8184e5d2b2eb21b5df74471f68be` |

Export reads committed review revision `55e207e33551475b70cb69d172d8b788429e0577`
and workflow revision `c98b11fb02843cdd6fd19af2a52ffb2e2661ee3a`. Main text and
render-time ledger use the exact `b74494376002dc87094aeb5ad733c9b0ea45cdd4`
snapshot, not current progress text. New CI findings are not substituted into
those historical render inputs. The exporter and its bundled guide are the
committed workflow snapshot; later documentation is not inserted retroactively.

All 67 tar members pass canonical/unique-path, regular-file/directory, mode and
size checks before data-filtered extraction to a fresh `/tmp` directory outside
the checkout. Isolated `/usr/bin/python3 -I -B` verification with PATH `/no-git`,
HOME `/nonexistent` and no scientific packages gives exactly the original result.
A second isolated verification installs a Python-event guard rejecting reads
under the original checkout and subprocess execution. Its explicit checkout-read
canary is rejected; verification then records zero original-path events. This is
a Python-level check, not OS containment, hermetic installation or another host.
The temporary extraction is removed; original archive and component remain.

After safe transfer/extraction, supply the independently retained index digest:

```sh
python3 -I -B /relocated/component/benchmark_tools/bundle_publication_review.py \
  verify /relocated/component \
  --manifest-sha256 d2b696bdeb9d70f13e1d5b7ff6e52808c0267f5f2f28b4e56049b654e2d5557a
```

Use the relative HTML entrypoint for portable direct-link navigation. The
byte-preserved PDF may retain workstation-specific link annotations. The receipt
retains build/verify argv and exact audit-guard source for repeatable readback.

## Remaining Publication Work

This closes the seven-page review-copy archive gap, not the complete executable
study release. Linked documents' transitive evidence, raw inputs and complete
runtime dependencies are excluded. Archive integrity does not clear data rights,
validate scientific claims, establish controlled timing or deposit a public release.

At 04:27:17 UTC the existing source-c98 CI has five live test jobs and successful
wheel/docs jobs. No new test log/outcome is inferred from setup or queued status;
no restart occurs. Other-QfO uncertainty, TreeFam sources, controlled-resource
evidence, transitive dependency/rights review, final manuscript reconciliation,
versioned release and public deposition remain open. Timing stays deferred
without contention polls, quiet-window questions or DGX/workload/service actions.
