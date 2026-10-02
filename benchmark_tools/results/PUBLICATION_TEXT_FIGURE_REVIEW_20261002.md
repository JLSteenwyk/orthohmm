# Main Text With Figure Appendix

The [accepted working-review PDF](publication_main_with_figures_20261002_v3/document.pdf)
contains the original nine-page manuscript, four caption/guide pages and
14 existing vector figure pages. Readers can inspect text and graphics without
opening separate figure files. Six original figure links now navigate internally;
the guide has 14 links and the PDF has 16 bookmarks.

The appendix covers the method, current OrthoBench/Three Kingdoms overview,
individual corrected QfO endpoints, paired OrthoBench uncertainty, the corrected
SwissTrees factorial, GO/EC eligible-pair decomposition, matched-recall search,
synthetic null scores, fixed/variable-length simulations, species-tree controls,
the complete parameter neighborhood, YGOB overlap strata and the WGD application.
No controlled-resource figure is invented; timing is still deferred.

## Actual Verification

- Each figure PDF matches its producer's retained output hash. Main PDF and
  25 of 28 figure/manifest files match Git at source `a0ad83ab`.
- Three direct sources are explicitly local-only: the fixed-length PDF/manifest
  and the large parameter manifest. Their retained output/result or committed
  compact-readback bindings are checked. Available direct result-file hashes
  are checked, not all transitive scientific provenance or statistical arithmetic.
- All 23 original text/figure pages preserve exact dimensions, text-word geometry
  and rendered pixels. No plot, bootstrap, inference or scientific score is rerun.
- Separate full-document readback verifies all original manuscript links, six
  intended redirects, 14 guide links, 16 bookmarks, captions and 62 evidence pins.
  The saved PDF needs no repair and emits no MuPDF warning.
- Guide text has zero detected bounds violations or block overlaps.
- All 14 figures and four guide layouts were visually inspected in the first
  assembly. Final guide page 10 is reinspected after caption correction; the
  other appended layouts are reused through pixel identity. The original main
  visual review is reused, not newly repeated, through identical page pixels.
- A negative CLI attempt refuses an existing output directory; PDF, assembly
  report and source bytes remain unchanged.

The [assembly](publication_main_with_figures_20261002_v3/assembly.json),
[separate readback](publication_main_with_figures_20261002_v3/readback.json) and
[manual visual closure](publication_main_with_figures_20261002_v3/visual_review.json)
record different scopes. The first two deliberately do not claim automatic
visual review. Readback uses the same PyMuPDF renderer, not independent
cross-viewer or formal PDF certification. qpdf is not installed.

## Preserved Attempts

The initial source check fails before output creation because the fixed-length
figure is not in Git. Distinguish exactly three local-only sources rather than
silently calling them committed. The first assembly emits six transient xref
warnings during deletion/reinsertion of links; its saved PDF reopens cleanly.
Visual review also finds guide A2 incorrectly mentioning a QfO mean absent
from the overview. Preserve the first assembly as unaccepted.

V2 corrects that caption and avoids deletion, but its independent link-signature
check fails before readback output: page copying converted original file URI
actions into remote-PDF actions. A bounded in-memory check demonstrates restoring
the original URI. V3 restores all 45 such actions and redirects only included
figures; the unchanged strict comparison passes. Earlier sources/reports/PDFs
remain retained, not overwritten or promoted.

## Artifact Pins

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| PDF | 470,650 | `2d225baf0545bf95ad28ce56a3afe735cc27c1a86effc0d6a85fdc0d73d7efb1` |
| Assembly | 56,584 | `861b869a4e359106843692c50a3e43e7805c177b32cce151fcf9cc2ab271176e` |
| Readback | 19,927 | `3b3b90af277d50e85cd06e35f7c25c48d532daf55b1184349141c40ed762b182` |
| Visual closure | 8,184 | `a5107025c2c8e9cf2ff8ce351e9985e31185899528b536825f098383018aa112` |

The [assembly source](publication_figure_appendix_20261002_v3.py) accepts a fresh
`--output` directory. It requires retained original PDFs/manifests and the
local predecessor receipt; its closed source is run provenance, not a portable
inference/release recipe. Rendered PNGs, large local manifests and rejected
prototypes are retained locally, not added as raw payloads to Git.

This is self-contained for the included text/figures, not for audit/data/code
links or executable reproduction. Original figure page sizes vary; final
journal typesetting, controlled resources, other-QfO uncertainty and complete
publication/release reconciliation remain open. No scientific default, score,
timing admission, source/runtime recipe, scientific or timing-helper implementation, scheduler,
service or unrelated-job state changes. No DGX, contention poll, archive rebuild
or new inference. The original publication goal remains incomplete.
